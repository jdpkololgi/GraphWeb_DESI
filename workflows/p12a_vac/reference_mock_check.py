"""Separate full-BRIGHT version check from BAO-selected reference-mock counts."""
import argparse,json,os,hashlib
from pathlib import Path
import numpy as np
import fitsio,h5py,hdf5plugin,healpy as hp
from mock_selection_test import cap
ROOT=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2')
REPO=Path(__file__).resolve().parents[2]
ZE=np.linspace(.1,.6,51)

def chunks(p,cols):
 if p.suffix=='.h5':
  with h5py.File(p) as f:
   t=f['LSS']
   for s in range(0,len(t[cols[0]]),250000):yield {k:t[k][s:s+250000] for k in cols}
 else:
  with fitsio.FITS(p) as f:
   n=f[1].get_nrows()
   for s in range(0,n,250000):yield f[1].read(rows=np.arange(s,min(s+250000,n)),columns=cols)

def hist(p,common,quality=False,weights=False):
 cols=['RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE'] if quality else ['RA','DEC','Z']
 if weights:cols+=['WEIGHT']
 n=np.zeros((2,50),dtype='i8');w=np.zeros((2,50));nall=0
 for d in chunks(p,cols):
  z=d['Z_not4clus' if quality else 'Z'];k=common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
  if quality:k&=(d['ZWARN']==0)&(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
  c=cap(d['RA'],d['DEC']);nall+=int(k.sum())
  for i in (0,1):
   keep=k&(c==i);n[i]+=np.histogram(z[keep],ZE)[0]
   w[i]+=np.histogram(z[keep],ZE,weights=d['WEIGHT'][keep] if weights else None)[0]
 print(p.name,'selected',nall,flush=True)
 return {'count':n.tolist(),'weighted':w.tolist(),'source':str(p),'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}

def main(out):
 out.mkdir(parents=True,exist_ok=True)
 if (out/'COMPLETE.json').exists():raise FileExistsError('completed run exists')
 common=np.load('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_mock_test_ph006_20260926_v1/common.npy')
 data={}
 p=ROOT/'LSS/loa-v1/LSScats/v1.1'
 data['loa_v1p1_full_bright']=hist(p/'BGS_BRIGHT_full_HPmapcut.dat.fits',common,quality=True)
 for c in ['SGC','NGC']:data['loa_v1p1_bao_'+c]=hist(p/f'nonKP/BGS_BRIGHT-21.35_{c}_clustering.dat.fits',common,weights=True)
 for phase in range(2,7):
  base=ROOT/f'mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl{phase}'
  for version,tracer,ext in [('kibo-v1','BGS_ANY-02','fits'),('loa-v1','BGS_BRIGHT-02','h5')]:
   for c in ['SGC','NGC']:
    path=base/f'{version}/mock{phase}/LSScats/{tracer}_{c}_clustering.dat.{ext}'
    key=f'ph{phase:03d}_{version}_{tracer}_{c}'
    data[key]=hist(path,common,weights=True) if path.exists() else {'source':str(path),'missing':True}
 # Does the Loa-named -02 product itself contain faint targets?
 faint={}
 for version,tracer,ext in [('kibo-v1','BGS_ANY-02','fits'),('loa-v1','BGS_BRIGHT-02','h5')]:
  path=ROOT/f'mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl6/{version}/mock6/LSScats/{tracer}_full_HPmapcut.dat.{ext}'
  cnt={};shell=np.zeros((2,2,5),dtype='i8')
  for d in chunks(path,['RA','DEC','Z_not4clus','ZWARN','BGS_TARGET']):
   k=(d['ZWARN']==0)&common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)];c=cap(d['RA'],d['DEC']);b=(d['BGS_TARGET']&2)!=0
   for flag,nn in zip(*np.unique(d['BGS_TARGET'][k],return_counts=True)):cnt[str(int(flag))]=cnt.get(str(int(flag)),0)+int(nn)
   for i in (0,1):
    for j in (0,1):shell[i,j]+=np.histogram(d['Z_not4clus'][k&(c==i)&(b==bool(j))],np.linspace(.1,.6,6))[0]
  faint[version]={'source':str(path),'successful_common_flags':cnt,'counts_cap_nonbright_bright_shell':shell.tolist()}
 result={'z_edges':ZE.tolist(),'catalogues':data,'faint_target_audit':faint,'job':os.environ.get('SLURM_JOB_ID'),'phases':[2,3,4,5,6], 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'notes':['Only already exposed phases; no truth or reserved phases opened.','Full-BRIGHT v1.1 uses our frozen quality25+GALAXY rule; BAO products use their stored selection.','Fixed previous common occupied-pixel sky; no n(z) rescaling.','WEIGHT sums are reported separately, not interchangeable with unweighted counts or field inputs.']}
 (out/'COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n')
 print('COMPLETE',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.out)
