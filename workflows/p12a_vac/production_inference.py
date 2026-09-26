"""Restartable halo48 inference: preserve the golden-tested arithmetic and solver."""
from loa_trial import *

def verify_core(marker,ready):
    receipt=read(marker);root=marker.parent.parent
    if receipt['input_sha256']!=digest(root/'INPUTS_READY.json'):raise ValueError('resume input changed')
    expected={c.get('core_id',i):c for i,c in enumerate(ready['cases'])}
    if receipt['case']!=expected[receipt['core_id']]:raise ValueError('resume ownership changed')
    for key in ['fits','draws']:
        if digest(root/receipt[key])!=receipt[key+'_sha256']:raise ValueError('corrupt completed '+key)
    return receipt

def verify_complete(root):
    receipt=read(root/'VAC_COMPLETE.json')
    if digest(root/'INPUTS_READY.json')!=receipt['inputs_marker_sha256']:raise ValueError('complete input changed')
    if digest(receipt['catalogue'])!=receipt['catalogue_sha256']:raise ValueError('corrupt complete catalogue')
    for item in receipt['draw_shards']:
        if digest(item['path'])!=item['sha256']:raise ValueError('corrupt complete draws')
    return receipt

def infer(root):
    parity=read(root/'RESPONSE_PARITY.json')
    if not parity['distance_exact'] or not parity['support_exact']:raise ValueError('response parity gate failed')
    ready=read(root/'INPUTS_READY.json')
    for name,sha in ready['arrays_sha256'].items():
        if digest(root/name)!=sha:raise ValueError('input hash mismatch: '+name)
    if (root/'VAC_COMPLETE.json').exists():
        verify_complete(root);print('REUSED COMPLETE',root,flush=True);return
    selection=read(ready['selection_path']);p3=read(B/'ph006/p3_fields/field_manifest.json');schema=read(p3['frozen_schema']);spline=read(p3['ntilde_spline'])
    for grid in ready['grids'].values():
        if grid['padding_mpc']!=schema['grid']['padding_mpc']:raise ValueError('grid padding differs from frozen schema')
    catalogue=np.load(root/'catalogue.npz');xyz=catalogue['xyz'];cap=catalogue['cap'];shell=catalogue['shell'];z=catalogue['z']
    angular=np.load(root/'galaxy_angular_support.npy');random=np.load(root/'random_angular.npz');angles=np.load(root/'boundary_angles.npz')
    sr=read(C/'summaries/ph006/OOF_SUMMARY_COMPLETE.json');ck=torch.load(sr['checkpoint'],map_location='cuda',weights_only=False)
    if digest(sr['checkpoint'])!=sr['checkpoint_sha256']:raise ValueError('encoder changed')
    torch.set_num_threads(8);torch.backends.cudnn.allow_tf32=True;torch.backends.cuda.matmul.allow_tf32=False
    if digest(C/'posterior/fmpe_estimator.pt')!=read(C/'posterior/P12A_COMPLETE.json')['checkpoint_sha256']:raise ValueError('posterior checkpoint changed')
    model=u.UPatch().cuda().eval();model.load_state_dict(ck['state_dict']);posterior,pck=reconstruct_fmpe(C/'posterior/fmpe_estimator.pt','cuda')
    train=read(C/'dataset/P12A_DATASET_READY.json')
    with np.load(train['training']['path']) as f:lo=f['context'].min(axis=0);hi=f['context'].max(axis=0)
    outputs=[];draw_files=[];qa=[];seen=set()
    (root/'shards').mkdir(exist_ok=True)
    for local_index,case in enumerate(ready['cases']):
        k=case.get('core_id',local_index)
        marker=root/f'shards/core_{k:06d}_COMPLETE.json'
        if marker.exists():
            receipt=verify_core(marker,ready)
            outputs.append(fitsio.read(root/receipt['fits']))
            draw_files.append(dict(path=str(root/receipt['draws']),sha256=receipt['draws_sha256']))
            qa.append(receipt['qa']);seen.update(outputs[-1]['TARGETID'].tolist());continue
        name='NGC' if case['cap'] else 'SGC';grid=ready['grids'][name];origin=np.asarray(grid['origin_mpc']);start=np.asarray(case['core_start']);stop=np.asarray(case['core_stop']);shape=np.asarray(grid['shape'])
        frac=(xyz-origin)/5-.5;voxel=np.floor((xyz-origin)/5).astype('i8')
        rows=np.flatnonzero((cap==case['cap'])&(shell>=0)&np.all((voxel>=start)&(voxel<stop),axis=1))
        ids=catalogue['targetid'][rows]
        if seen.intersection(ids.tolist()):raise ValueError('output core ownership collision')
        seen.update(ids.tolist());cs=np.maximum(((start-48)//8)*8,0);ce=np.minimum(((stop+48+7)//8)*8,shape)
        context=np.flatnonzero((cap==case['cap'])&np.all((frac>=cs-1)&(frac<ce),axis=1))
        fields=obs.rebuild_fields(xyz[context],shell[context],grid,cs,ce,angular,256,schema,spline,selection,0,name)
        if not all(np.isfinite(fields[n]).all() for n in u.CHANNELS):raise ValueError('nonfinite input fields')
        patch=FieldPatch(k,-1,case['cap'],tuple(u.CHANNELS),np.stack([fields[n] for n in u.CHANNELS]),cs,ce,start,stop,tuple(slice(int(a),int(b)) for a,b in zip(start-cs,stop-cs)),rows,frac[rows],frac[rows]-cs,start-cs,ce-stop)
        base=predict(model,patch,ck)
        distance,support=response_at(xyz[rows],grid,random['support']&((random['domain']//2)==case['cap']),angles[name],selection)
        nt=ntilde_at_rows(selection,cap[rows],z[rows].astype('f4'))
        x=np.column_stack([base,z[rows].astype('f4'),np.log(nt),cap[rows],np.log1p(distance)]).astype('f4')
        outside=np.any((x<lo)|(x>hi),axis=1);eligible=support&np.isfinite(x).all(axis=1)
        draws=np.full((len(rows),512,3),np.nan,dtype='f4')
        torch.manual_seed(20260925+k);torch.cuda.manual_seed_all(20260925+k)
        if eligible.any():
            scaled=((x[eligible]-pck['context_mean'])/pck['context_std']).astype('f4');sample=sample_posterior(posterior,scaled,512,128,'cuda')
            draws[eligible]=theta_to_eigenvalues(sample*np.asarray(pck['theta_std'])+np.asarray(pck['theta_mean'])).astype('f4')
            if not np.isfinite(draws[eligible]).all() or np.any(np.diff(draws[eligible],axis=-1)<0):raise ValueError('invalid posterior draws')
        dtype=[('TARGETID','i8'),('RA','f8'),('DEC','f8'),('Z','f8'),('CAP','i1'),('CORE_ID','i4'),('QUALITY','u2'),('SUPPORTED','?'),('BASE_EIGENVALUES','f4',(3,)),('EIGENVALUE_MEAN','f4',(3,)),('EIGENVALUE_Q05','f4',(3,)),('EIGENVALUE_Q16','f4',(3,)),('EIGENVALUE_Q50','f4',(3,)),('EIGENVALUE_Q84','f4',(3,)),('EIGENVALUE_Q95','f4',(3,)),('P_WEB','f4',(4,)),('BOUNDARY_MPC','f4'),('NTILDE_MPC3','f4')]
        out=np.zeros(len(rows),dtype=dtype)
        for n,v in [('TARGETID',ids),('RA',catalogue['ra'][rows]),('DEC',catalogue['dec'][rows]),('Z',z[rows]),('CAP',cap[rows]),('CORE_ID',k),('SUPPORTED',eligible),('BASE_EIGENVALUES',base),('BOUNDARY_MPC',distance),('NTILDE_MPC3',nt)]:out[n]=v
        out['QUALITY']=np.uint16(32)+(~eligible).astype('u2')+2*(shell[rows]==3).astype('u2')+4*(distance<10.345846881466155).astype('u2')+8*(distance<20.69169376293231).astype('u2')+16*outside.astype('u2')
        for n in ['EIGENVALUE_MEAN','EIGENVALUE_Q05','EIGENVALUE_Q16','EIGENVALUE_Q50','EIGENVALUE_Q84','EIGENVALUE_Q95','P_WEB']:out[n]=np.nan
        if eligible.any():
            out['EIGENVALUE_MEAN'][eligible]=draws[eligible].mean(axis=1);q=np.quantile(draws[eligible],[.05,.16,.5,.84,.95],axis=1)
            for j,n in enumerate(['EIGENVALUE_Q05','EIGENVALUE_Q16','EIGENVALUE_Q50','EIGENVALUE_Q84','EIGENVALUE_Q95']):out[n][eligible]=q[j]
            classes=np.sum(draws[eligible]>.2,axis=-1)
            out['P_WEB'][eligible]=np.stack([(classes==j).mean(axis=1) for j in range(4)],axis=1)
            if not np.allclose(out['P_WEB'][eligible].sum(axis=1),1):raise ValueError('class probabilities not normalized')
        file=root/f'shards/core_{k:06d}_draws.npz';tmp=file.with_suffix('.tmp.npz');np.savez_compressed(tmp,TARGETID=ids,eigenvalue_draws=draws);tmp.replace(file);draw_files.append(dict(path=str(file),sha256=digest(file)))
        outputs.append(out);qa.append(dict(**case,rows=len(rows),supported=int(eligible.sum()),outside_training_envelope=int(outside.sum()),context_rows=len(context),raw_input_count_sum=float(fields['counts'].sum(dtype='f8')),median_base= np.median(base,axis=0).tolist()))
        core_fits=root/f'shards/core_{k:06d}.fits';tmp=core_fits.with_suffix('.tmp.fits');fitsio.write(tmp,out,extname='ENV_POSTERIOR',clobber=True);tmp.replace(core_fits)
        save(marker,dict(core_id=k,case=case,input_sha256=digest(root/'INPUTS_READY.json'),fits=str(core_fits.relative_to(root)),fits_sha256=digest(core_fits),draws=str(file.relative_to(root)),draws_sha256=digest(file),qa=qa[-1]))
        print('core',k,'rows',len(rows),'supported',eligible.sum(),'outside',outside.sum(),flush=True)
    out=np.concatenate(outputs);out.sort(order='TARGETID');file=root/'VAC_SHARD.fits'
    tmp=file.with_suffix('.tmp.fits')
    fitsio.write(tmp,out,extname='ENV_POSTERIOR',header={'PROVIS':True,'MODEL':'P12A_H48','NSAMPLE':512,'LTHRESH':.2,'RSMOOTH':7.,'RUNIT':'Mpc/h','ZTARGET':.2},clobber=True)
    tmp.replace(file)
    reread=fitsio.read(file);assert np.array_equal(reread['TARGETID'],out['TARGETID'])
    save(root/'VAC_COMPLETE.json',dict(schema='desi-loa-p12a-halo48-canary-vac-v1',technical_complete=True,science_release_ready=False,scope=ready.get('scope','production benchmark'),rows=len(out),supported_rows=int(out['SUPPORTED'].sum()),catalogue=str(file),catalogue_sha256=digest(file),draw_shards=draw_files,cases=qa,inputs_marker_sha256=digest(root/'INPUTS_READY.json'),posterior_checkpoint_sha256=digest(C/'posterior/fmpe_estimator.pt'),quality_bits={'1':'unsupported; posterior null','2':'sparse shell z>=0.45','4':'boundary distance <R/h','8':'boundary distance <2R/h','16':'one or more posterior features outside training min/max envelope','32':'provisional diagnostic only'},class_order=['void','sheet','filament','knot'],job=os.environ['SLURM_JOB_ID']))
    print('VAC COMPLETE',len(out),flush=True)

