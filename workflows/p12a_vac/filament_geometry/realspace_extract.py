"""Bounded three-region extraction for a native-matter / RSD / sampling audit."""
import json,os,argparse
from pathlib import Path
import numpy as np
import fitsio
from astropy.cosmology import Planck18
from prepare_parent import digest,coordinates
from mock_truth import ROOT

RAW=Path('/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v0.1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph006.fits')
PARENT=ROOT/'ph006/catalogues/annotated_parent/ph006_bgs_bright_parent_with_tweb_eigs_rs7_ngrid2048_thr0p2_15d.fits'

def keys(d):
 a=np.empty(len(d),dtype=[('RA','f8'),('DEC','f8'),('Z','f4')])
 for k in a.dtype.names:a[k]=d[k]
 return a

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 assert os.environ.get('SLURM_JOB_ID')
 if a.out.exists():raise FileExistsError(a.out)
 examples=Path(__file__).resolve().parents[3]/'docs/evidence/filament_mock_examples_20261017/EXAMPLES.json';m=json.loads(examples.read_text());targets=[v['targetid'] for v in m['examples']]
 with fitsio.FITS(PARENT) as f:
  centers=np.concatenate([f[1].read(rows=[t-1]) for t in targets])
 np.testing.assert_array_equal(centers['TARGETID'],targets)
 zgrid=np.linspace(0,1,20001);rgrid=Planck18.comoving_distance(zgrid).value
 def xyz(d,z):return coordinates(d['RA'],d['DEC'],np.interp(d[z],zgrid,rgrid))
 center_real=xyz(centers,'Z_COSMO');center_rs=xyz(centers,'Z')
 parts={name:[[] for _ in targets] for name in ['raw','parent']};sources=[]
 for name,path in [('parent',PARENT),('raw',RAW)]:
  with fitsio.FITS(path) as f:
   total=f[1].get_nrows()
   for lo in range(0,total,500000):
    d=f[1][lo:min(lo+500000,total)]
    real=xyz(d,'Z_COSMO');rs=xyz(d,'Z')
    for i in range(len(targets)):
     keep=(np.sum((real-center_real[i])**2,axis=1)<45**2)|(np.sum((rs-center_rs[i])**2,axis=1)<45**2)
     if keep.any():parts[name][i].append(d[keep])
    if lo%10000000==0:print(name,lo,total,flush=True)
   sources.append({'path':str(path),'rows':total,'bytes':path.stat().st_size,'mtime_ns':path.stat().st_mtime_ns,'full_hash_checked':False})
 a.out.mkdir(parents=True)
 ci=np.load(ROOT/'ph006/p1_canonical/canonical_index.npz');observed=np.sort(ci['targetid'])
 from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
 manifest={'job':os.environ['SLURM_JOB_ID'],'sources':sources,'examples_sha256':digest(examples),'script_sha256':digest(__file__),'regions':[]}
 for i,row in enumerate(m['examples']):
  raw=np.concatenate(parts['raw'][i]);parent=np.concatenate(parts['parent'][i]);pk=keys(parent);rk=keys(raw);order=np.argsort(pk);at=np.searchsorted(pk[order],rk);match=(at<len(pk))&(pk[order[np.minimum(at,len(pk)-1)]]==rk)
  ids=np.full(len(raw),-1,dtype='i8');ids[match]=parent['TARGETID'][order[at[match]]]
  ic=np.flatnonzero(ids==targets[i]);assert len(ic)==1
  cen=int(ic[0]);real=xyz(raw,'Z_COSMO');rs=xyz(raw,'Z');native=np.column_stack([raw[k] for k in ['X_COM','Y_COM','Z_COM']]);rep=np.column_stack([raw[k] for k in ['X_REP','Y_REP','Z_REP']])
  assert len(np.unique(ids[match]))==match.sum()
  idx=np.searchsorted(observed,ids);obs=(ids>0)&(idx<len(observed))&(observed[np.minimum(idx,len(observed)-1)]==ids)
  halo_path=Path('/global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/AbacusSummit_base_c000_ph006/halos/z0.200/halo_info')/f"halo_info_{int(raw['FILE_NUM'][cen]):03d}.asdf"
  cat=CompaSOHaloCatalog(halo_path,fields=['x_com'],subsamples=False,cleaned=False,convert_units=True,verbose=False);host=np.asarray(cat.halos['x_com'][int(raw['HALO_INDEX'][cen])]);del cat
  # Geometry-only local transformation diagnostic. No labels used to fit it.
  dn=(native-native[cen]+1000)%2000-1000;dr=real-real[cen]
  nearby=(np.linalg.norm(dr,axis=1)<30)&(np.linalg.norm(dn,axis=1)<30)&(np.linalg.norm(dr,axis=1)>1)
  train=np.flatnonzero(nearby)[::2];test=np.flatnonzero(nearby)[1::2]
  transform=np.linalg.lstsq(dn[train],dr[train],rcond=None)[0];resid=np.linalg.norm(dn[test]@transform-dr[test],axis=1)
  result=dict(row,key=row['key'],raw_rows=len(raw),parent_rows=len(parent),parent_matches=int(match.sum()),central_raw_index=cen,host_native=host.tolist(),galaxy_native=native[cen].tolist(),host_offset_mpc_h=float(np.linalg.norm((host-native[cen]+1000)%2000-1000)),local_transform_native_to_observer=transform.tolist(),transform_singular_values=np.linalg.svd(transform,compute_uv=False).tolist(),transform_test_rows=len(test),transform_residual_Mpc_quantiles=np.quantile(resid,[0,.5,.95,1]).tolist(),central_source={k:float(raw[k][cen]) for k in ['RA','DEC','Z','Z_COSMO','R_MAG_APP','FILE_NUM','HALO_INDEX','BOX_INDEX']},native_bounds=[native.min(0).tolist(),native.max(0).tolist()],halo_path=str(halo_path))
  path=a.out/(row['key']+'.npz');np.savez_compressed(path,raw=raw,parent=parent,targetid=ids,observed=obs,real_Mpc=real,rs_Mpc=rs,native_mpc_h=native,rep=rep,central_index=cen,host_native=host,transform=transform)
  result['sha256']=digest(path);manifest['regions'].append(result);print(json.dumps(result),flush=True)
 (a.out/'EXTRACTION.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
