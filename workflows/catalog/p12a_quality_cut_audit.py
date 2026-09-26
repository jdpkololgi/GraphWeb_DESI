"""Census existing DESI quality flags, without altering catalogue selection."""
import argparse, hashlib, json, os, socket
from pathlib import Path
import fitsio
import numpy as np


def masks(d, zcol):
    good = d['ZWARN'] == 0
    spec = np.char.strip(d['SPECTYPE'].astype('U')) == 'GALAXY'
    finite = np.isfinite(d[zcol]) & (d[zcol] > 0)
    return dict(zwarn0=good, legacy25=good & (d['DELTACHI2'] >= 25) & spec,
                lss40=good & (d['DELTACHI2'] > 40), galaxy=spec,
                delta25_to40=good & (d['DELTACHI2'] >= 25) & (d['DELTACHI2'] <= 40),
                delta_equal40=good & (d['DELTACHI2'] == 40), finite_positive_z=finite)


def census(path, zcol):
    before=path.stat(); totals={}; dig=hashlib.sha256()
    with fitsio.FITS(path) as f:
        n=f[1].get_nrows()
        for start in range(0,n,500000):
            d=f[1].read(rows=range(start,min(start+500000,n)),columns=['TARGETID',zcol,'ZWARN','DELTACHI2','SPECTYPE'])
            dig.update(d.tobytes()); m=masks(d,zcol)
            m['legacy_only']=m['legacy25'] & ~m['lss40'];m['lss_only']=m['lss40'] & ~m['legacy25']
            for name,lo,hi in [('all',-np.inf,np.inf),('active',.15,.55),('shell0',.15,.25),('shell1',.25,.35),('shell2',.35,.45),('shell3',.45,.55)]:
                region=np.ones(len(d),bool) if name=='all' else ((d[zcol]>=lo)&(d[zcol]<hi))
                out=totals.setdefault(name,dict(rows=0,**{k:0 for k in m}))
                out['rows']+=int(region.sum())
                for k,v in m.items():out[k]+=int((v & region).sum())
    after=path.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValueError('source changed')
    return dict(path=str(path),bytes=before.st_size,mtime_ns=before.st_mtime_ns,scanned_columns_sha256=dig.hexdigest(),counts=totals)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):raise RuntimeError('compute allocation required')
    sources=[('existing_graphweb',Path('/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits'),'Z'),('loa_full',Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits'),'Z_not4clus')]
    records={}
    for label,path,zcol in sources:
        records[label]=census(path,zcol);print(label,records[label]['counts'],flush=True)
    r=dict(schema='p12a-quality-cut-census-v1',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),job=os.environ['SLURM_JOB_ID'],catalogues=records,selection_adopted=False,ready_for_desi_canary=False)
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
