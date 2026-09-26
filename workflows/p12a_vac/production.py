"""Frozen full-survey planning, independently restartable shards, and checked merge."""
import argparse,fcntl,shutil,subprocess
from production_inference import *

INPUT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')

def bind_inputs(dest,ready):
    dest.mkdir(parents=True,exist_ok=True)
    for name in [*ready['arrays_sha256'],'RESPONSE_PARITY.json']:
        if not (dest/name).exists():(dest/name).symlink_to(INPUT/name)
    save(dest/'INPUTS_READY.json',ready)

def plan(root):
    if (root/'PLAN.json').exists():raise FileExistsError('plan exists')
    ready=read(INPUT/'INPUTS_READY.json')
    for name,value in ready['arrays_sha256'].items():
        if digest(INPUT/name)!=value:raise ValueError('input changed')
    d=np.load(INPUT/'catalogue.npz');xyz=d['xyz'];cap=d['cap'];shell=d['shell'];cases=[]
    for c,name in [(0,'SGC'),(1,'NGC')]:
        grid=ready['grids'][name];rows=np.flatnonzero((cap==c)&(shell>=0));voxel=np.floor((xyz[rows]-grid['origin_mpc'])/5).astype('i8');cores,counts=np.unique(voxel//20,axis=0,return_counts=True)
        for core,n in zip(cores,counts):
            cases.append(dict(core_id=len(cases),cap=c,core_start=(core*20).tolist(),core_stop=np.minimum(core*20+20,grid['shape']).tolist(),active_rows=int(n)))
    if sum(c['active_rows'] for c in cases)!=ready['active_rows']:raise ValueError('ownership census')
    groups=[[] for _ in range(128)];cost=np.zeros(128)
    for case in sorted(cases,key=lambda c:-(c['active_rows']+200)):
        j=int(np.argmin(cost));groups[j].append(case);cost[j]+=case['active_rows']+200
    root.mkdir(parents=True,exist_ok=True);(root/'logs').mkdir(exist_ok=True)
    for j,group in enumerate(groups):
        item=dict(ready,cases=sorted(group,key=lambda c:c['core_id']),scope='full-survey production shard; provisional science status')
        bind_inputs(root/f'parts/{j:03d}',item)
    # Truth-free dense benchmark plus the previous canary regions.
    keys={(c['cap'],tuple(c['core_start'])) for c in ready['cases']}
    bench=[c for c in cases if (c['cap'],tuple(c['core_start'])) in keys]
    for cap_id in [0,1]:
        dense=max([c for c in cases if c['cap']==cap_id],key=lambda c:c['active_rows'])
        if dense not in bench:bench.append(dense)
    bind_inputs(root/'benchmark',dict(ready,cases=bench,scope='dense and prior canary benchmark'))
    save(root/'PLAN.json',dict(schema='p12a-loa-full-survey-v1',science_rows=ready['active_rows'],context_rows=ready['context_rows'],cores=len(cases),parts=len(groups),part_input_sha256={str(j):digest(root/f'parts/{j:03d}/INPUTS_READY.json') for j in range(len(groups))},input_root=str(INPUT),input_sha256=digest(INPUT/'INPUTS_READY.json'),draws=512,seed='20260925 + stable global core_id',core_voxels=20,halo_voxels=48,alignment=8,benchmark_rows=sum(c['active_rows'] for c in bench),benchmark_cores=len(bench),science_release_ready=False,selection_shift='Count-density excess increases with redshift, approximately 25% at z=.35-.45 and 90% at .45-.55. Frozen selection retained; systematic uncertainty not marginalized.'))
    print(json.dumps(read(root/'PLAN.json'),indent=2),flush=True)

def check_sources(root):
    frozen=read(root/'SOURCE_MANIFEST.json')
    for name,value in frozen['files'].items():
        if digest(root/'source'/name)!=value:raise ValueError('frozen source changed: '+name)
    for path,value in frozen['artifacts'].items():
        if digest(path)!=value:raise ValueError('frozen artifact changed: '+path)
    if digest(root/'PLAN.json')!=frozen['plan_sha256']:raise ValueError('plan changed')
    science=Path(frozen['science_source']);manifest=read(science.parent/'RUN_MANIFEST.json')
    for name,value in manifest['source_sha256'].items():
        if digest(science/name)!=value:raise ValueError('frozen science source changed: '+name)

def worker(root,index):
    check_sources(root);p=read(root/'PLAN.json');part=root/f'parts/{index:03d}'
    if digest(part/'INPUTS_READY.json')!=p['part_input_sha256'][str(index)]:raise ValueError('part changed')
    with (part/'worker.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        infer(part)


def validate_rows(rows,draws):
    yes=rows['SUPPORTED'];d=draws[yes]
    if not np.isfinite(d).all() or np.any(np.diff(d,axis=-1)<0):raise ValueError('invalid draws')
    if not np.isnan(draws[~yes]).all():raise ValueError('unsupported draws nonnull')
    if not np.array_equal((rows['QUALITY']&1)==0,yes):raise ValueError('support flag')
    if not np.all((rows['QUALITY']&32)!=0):raise ValueError('missing provisional flag')
    names=['EIGENVALUE_MEAN','EIGENVALUE_Q05','EIGENVALUE_Q16','EIGENVALUE_Q50','EIGENVALUE_Q84','EIGENVALUE_Q95','P_WEB']
    if any(not np.isnan(rows[name][~yes]).all() for name in names):raise ValueError('unsupported summary nonnull')
    if yes.any():
        vals=[d.mean(axis=1),*np.quantile(d,[.05,.16,.5,.84,.95],axis=1).astype('f4'),np.stack([((d>.2).sum(axis=-1)==j).mean(axis=1) for j in range(4)],axis=1).astype('f4')]
        for name,value in zip(names,vals):
            if not np.array_equal(rows[name][yes],value):raise ValueError('saved draw summary mismatch '+name)
        if not np.allclose(rows['P_WEB'][yes].sum(axis=1),1):raise ValueError('probability normalization')

def merge(root):
    check_sources(root);p=read(root/'PLAN.json')
    if (root/'FULL_VAC_COMPLETE.json').exists():
        c=read(root/'FULL_VAC_COMPLETE.json')
        if digest(c['catalogue'])!=c['catalogue_sha256']:raise ValueError('corrupt merged catalogue')
        print('FULL VAC already complete');return
    tables=[];draw_index=[]
    for j in range(p['parts']):
        part=root/f'parts/{j:03d}';ready=read(part/'INPUTS_READY.json')
        if digest(part/'INPUTS_READY.json')!=p['part_input_sha256'][str(j)]:raise ValueError('part input changed')
        complete=verify_complete(part);table=fitsio.read(complete['catalogue']);table.sort(order='TARGETID');visited=[]
        for case in ready['cases']:
            cid=case['core_id'];record=verify_core(part/f'shards/core_{cid:06d}_COMPLETE.json',ready)
            with np.load(part/record['draws']) as f:
                ids=f['TARGETID'];ix=np.searchsorted(table['TARGETID'],ids);rows=table[ix]
                if not np.array_equal(rows['TARGETID'],ids) or not np.all(rows['CORE_ID']==cid):raise ValueError('draw identity or owner')
                validate_rows(rows,f['eigenvalue_draws']);visited.extend(ids.tolist())
            draw_index.append(dict(core_id=cid,path=str(part/record['draws']),sha256=record['draws_sha256'],rows=len(ids)))
        if not np.array_equal(np.sort(visited),table['TARGETID']):raise ValueError('shard census')
        tables.append(table);print('validated part',j,flush=True)
    out=np.concatenate(tables);del tables;out.sort(order='TARGETID');d=np.load(INPUT/'catalogue.npz');expected=np.sort(d['targetid'][d['shell']>=0])
    if not np.array_equal(out['TARGETID'],expected) or not np.all(np.diff(expected)>0):raise ValueError('full TARGETID census failed')
    out['QUALITY'] |= np.where(out['Z']>=.35,64,0).astype('u2')
    path=root/'DESI_LOA_P12A_HALO48_FULL_SURVEY_VAC.fits';tmp=path.with_suffix('.tmp.fits')
    fitsio.write(tmp,out,extname='ENV_POSTERIOR',header={'PROVIS':True,'MODEL':'P12A_H48','NSAMPLE':512,'LTHRESH':.2,'RSMOOTH':7.,'RUNIT':'Mpc/h','ZTARGET':.2,'SELSHIFT':True},clobber=True);tmp.replace(path)
    if not np.array_equal(fitsio.read(path,columns=['TARGETID'])['TARGETID'],expected):raise ValueError('FITS readback')
    save(root/'DRAW_INDEX.json',dict(shards=draw_index));save(root/'FULL_VAC_COMPLETE.json',dict(schema=p['schema'],technical_complete=True,science_release_ready=False,rows=len(out),supported_rows=int(out['SUPPORTED'].sum()),quality_bit_counts={str(b):int(((out['QUALITY']&b)!=0).sum()) for b in [1,2,4,8,16,32,64]},quality_bits={'1':'unsupported/null','2':'sparse z>=.45','4':'boundary<R/h','8':'boundary<2R/h','16':'outside training feature envelope','32':'provisional; science release not qualified','64':'unresolved selection mismatch in z>=.35 shells; no correction or uncertainty marginalization'},catalogue=str(path),catalogue_sha256=digest(path),draw_index_sha256=digest(root/'DRAW_INDEX.json'),plan_sha256=digest(root/'PLAN.json'),source_manifest_sha256=digest(root/'SOURCE_MANIFEST.json'),selection_shift=p['selection_shift'],unique_complete_targetid_census=True,saved_draw_summary_validation=True))
    print('FULL VAC COMPLETE',len(out),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['plan','benchmark','worker','merge']);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--part',type=int);a=parser.parse_args()
    if a.stage=='plan':plan(a.root)
    elif a.stage=='benchmark':infer(a.root/'benchmark')
    elif a.stage=='worker':worker(a.root,a.part)
    else:merge(a.root)
