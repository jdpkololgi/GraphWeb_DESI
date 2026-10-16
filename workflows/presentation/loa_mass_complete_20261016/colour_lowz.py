"""Colour-audited empirical mass-complete Loa sample at 0.15 <= z < 0.20."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,socket
from pathlib import Path
import numpy as np
from run import OUT as SOURCE,PHOT,VAC,sha,write_json
from colour_split import masks
from colour_followup import analyse
OUT=SOURCE.parent/'loa_colour_lowz_20261016_v1'
def colour_limits():
    """Independent colour-population completeness audit from the Bright parent."""
    import fitsio
    b=np.load(SOURCE/'bright_parent.npz'); order=np.argsort(b['TARGETID']); ids=b['TARGETID'][order]
    col=np.full(len(ids),np.nan); seen=np.zeros(len(ids),bool)
    with fitsio.FITS(PHOT) as f:
        for start in range(0,f[1].get_nrows(),300000):
            a=f[1].read(rows=np.arange(start,min(start+300000,f[1].get_nrows())),columns=['TARGETID','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R','ZWARN','DELTACHI2','Z_not4clus'])
            at=np.searchsorted(ids,a['TARGETID']); ok=(at<len(ids))&(ids[np.minimum(at,len(ids)-1)]==a['TARGETID'])
            ok&=(a['ZWARN']==0)&(a['DELTACHI2']>40)&(a['Z_not4clus']>0)&(a['Z_not4clus']<.55)
            ix=order[at[ok]]; assert not seen[ix].any();seen[ix]=True
            x=a[ok]
            with np.errstate(divide='ignore',invalid='ignore'):
                c=-2.5*np.log10((x['FLUX_G']/x['MW_TRANSMISSION_G'])/(x['FLUX_R']/x['MW_TRANSMISSION_R']))
            valid=(x['FLUX_G']>0)&(x['FLUX_R']>0)&(x['MW_TRANSMISSION_G']>0)&(x['MW_TRANSMISSION_R']>0)
            col[ix]=np.where(valid,c,np.nan)
    assert seen.all()
    curve=json.loads((SOURCE/'mass_limits.json').read_text()); groups=masks(col); audit=[]
    offset=b['r']-b['rlim']
    for old in curve:
        if old['zlo']<.15 or old['zhi']>.20:continue
        rec=dict(old);rec['colour_populations']={};limits=[old['logmass_limit']]
        for group,mask in groups.items():
            ii=np.flatnonzero(mask&(b['z']>=old['zlo'])&(b['z']<old['zhi']))
            item={'rows':len(ii),'usable':len(ii)>=100}
            if item['usable']:
                faint=ii[offset[ii]>=np.quantile(offset[ii],.8)]
                values=b['mass'][faint]+.4*offset[faint]
                lim=float(np.quantile(values,.95)); limits.append(lim)
                # Independent sorted order-statistic reconstruction.
                v=np.sort(values); pos=.95*(len(v)-1); j=int(pos)
                check=v[j]+(pos-j)*(v[min(j+1,len(v)-1)]-v[j]);assert abs(check-lim)<1e-10
                item.update(p95=lim,faint_rows=len(faint))
            rec['colour_populations'][group]=item
        rec['usable']=old['usable'] and all(x['usable'] for x in rec['colour_populations'].values())
        rec['logmass_limit']=max(limits) if rec['usable'] else None;audit.append(rec)
    return audit

def main():
    import fitsio
    from numpy.lib.recfunctions import append_fields
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):raise RuntimeError('Compute allocation required')
    OUT.mkdir(exist_ok=False)
    curve=colour_limits();assert len(curve)==5 and all(b['usable'] for b in curve)
    limit=float(np.ceil(100*max(b['logmass_limit'] for b in curve))/100)
    write_json(OUT/'colour_mass_limits.json',curve)
    # The old variable-cut VAC is a lossless source for this stricter fixed cut:
    # every old bin limit in this window is <= the new common threshold.
    old=json.loads((SOURCE/'mass_limits.json').read_text())
    assert all(b['logmass_limit']<=limit for b in old if b['zlo']>=.15 and b['zhi']<=.20)
    a=fitsio.read(SOURCE/'LOA_P12A_MASS_COMPLETE_zdependent.fits')
    keep=(a['Z']>=.15)&(a['Z']<.20)&(np.log10(a['MASS_CG_15'])>=limit)
    a=a[keep];assert len(np.unique(a['TARGETID']))==len(a)
    groups=masks(a['GR_OBS']);counts={g:int(v.sum()) for g,v in groups.items()}
    hard={g:np.bincount(a['P_WEB'][v].argmax(1),minlength=4).tolist() for g,v in groups.items()}
    window={'zlo':.15,'zhi':.20,'logmass':limit,'mass_msun':float(10**limit),'rows':len(a),'counts':counts,'hard_counts':hard,'adequate_counts':min(counts.values())>=20000 and min(min(v) for v in hard.values())>=100}
    write_json(OUT/'windows.json',{'z015020':window});print(json.dumps(window),flush=True)
    results={'z015020':analyse(a)};write_json(OUT/'relationships.json',results)
    flag=np.full(len(a),-1,dtype='i2');flag[groups['blue']]=0;flag[groups['red']]=1
    sub=append_fields(a,['COLOUR_CLASS','COLOUR_LOGM_LIMIT'],[flag,np.full(len(a),limit)],usemask=False)
    fitsio.write(OUT/'LOA_COLOUR_MASS_COMPLETE_z015020.fits',sub,header={'MCFIX':limit,'MCZLO':.15,'MCZHI':.20,'COLRED':'0.75 < g-r < 1.5','COLBLUE':'g-r < 0.75'})
    ix=np.linspace(0,len(a)-1,257,dtype=int);fresh=fitsio.read(OUT/'LOA_COLOUR_MASS_COMPLETE_z015020.fits',rows=ix)
    original=fitsio.read(VAC,rows=a['SOURCE_ROW'][ix])
    for col in sub.dtype.names:
        eq=np.array_equal(fresh[col],sub[col][ix],equal_nan=True) if fresh[col].dtype.kind in 'fc' else np.array_equal(fresh[col],sub[col][ix]);assert eq,col
    for col in original.dtype.names:
        eq=np.array_equal(original[col],fresh[col],equal_nan=True) if original[col].dtype.kind in 'fc' else np.array_equal(original[col],fresh[col]);assert eq,col
    write_json(OUT/'MANIFEST.json',{'job':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'source_manifest_sha256':sha(SOURCE/'MANIFEST.json'),'code':{str(q):sha(q) for q in [Path(__file__),Path(__file__).with_name('run.py'),Path(__file__).with_name('colour_split.py'),Path(__file__).with_name('colour_followup.py')]},'outputs':{q.name:sha(q) for q in OUT.iterdir() if q.is_file()},'validation':'10 colour percentiles, every raw/controlled mean, 257-row FITS roundtrip and fresh all-original-column VAC parity passed; fixed threshold >= every old bin limit proves no source truncation','limitations':'Empirical flux completeness; observed colour; separate red/blue support; conditional sky bootstrap excludes SED/selection/model-transfer systematics.'})
if __name__=='__main__':main()
