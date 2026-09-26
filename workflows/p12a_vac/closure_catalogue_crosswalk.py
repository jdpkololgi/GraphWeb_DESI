"""Read-only Kibo/Loa data and ph006 mock product crosswalk."""
import json,os
from pathlib import Path
import numpy as np,fitsio,h5py,healpy as hp
import hdf5plugin  # Register compression filters used by official Loa mock HDF5.
import loa_trial as t
BASE=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2')
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_closure_20260926'
ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
def scan(path,real,common):
 counts=np.zeros((2,40),dtype='i8');full=np.zeros_like(counts);ntile=np.zeros((8,30),dtype='i8');ids=[];zs=[];tiles=set();rowtiles=set();stages={};n=0
 is_h5=path.suffix=='.h5';f=h5py.File(path) if is_h5 else fitsio.FITS(path);tab=f['LSS'] if is_h5 else f[1];keys=list(tab.keys()) if is_h5 else tab.get_colnames();nr=len(tab['TARGETID']) if is_h5 else tab.get_nrows()
 for start in range(0,nr,250000):
  sl=slice(start,min(start+250000,nr));d={k:tab[k][sl] for k in ['TARGETID','RA','DEC','Z_not4clus','ZWARN','NTILE','TILEID','TILES','GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED'] if k in keys} if is_h5 else tab[sl];n+=len(d['TARGETID']);z=d['Z_not4clus'];shell=np.searchsorted([.15,.25,.35,.45,.55],z,side='right')-1;valid=np.isfinite(z)&(shell>=0)&(shell<4);good=valid&(d['ZWARN']==0)
  if real:good&=(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
  cap=t.galactic_cap(d['RA'],d['DEC']);pix=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True);cm=common[pix];b=np.clip(cap*4+shell,0,7)
  for name,k in [('finite_science_z',valid),('zwarn0',valid&(d['ZWARN']==0)),('success',good),('common_success',good&cm)]:
   stages.setdefault(name,np.zeros(8,dtype='i8'));stages[name]+=np.bincount(b[k],minlength=8)
  for c in [0,1]:
   full[c]+=np.histogram(z[good&(cap==c)],np.linspace(.15,.55,41))[0];counts[c]+=np.histogram(z[good&cm&(cap==c)],np.linspace(.15,.55,41))[0]
  for j in range(8):ntile[j]+=np.bincount(np.minimum(d['NTILE'][good&cm&(b==j)],29).astype(int),minlength=30)
  ids.append(np.asarray(d['TARGETID'][good]));zs.append(np.asarray(z[good]));rowtiles.update(int(x) for x in np.unique(d['TILEID']))
  if 'TILES' in keys:
   for val in np.unique(d['TILES']):
    if isinstance(val,bytes):val=val.decode()
    for s in str(val).split('-'):
     if s.isdigit():tiles.add(int(s))
 f.close();ids=np.concatenate(ids);zs=np.concatenate(zs);order=np.argsort(ids);assert not np.any(np.diff(ids[order])==0)
 return dict(path=str(path),bytes=path.stat().st_size,rows=n,columns=keys,stages={k:v.tolist() for k,v in stages.items()},counts=counts.tolist(),fullcounts=full.tolist(),ntile_hist=ntile.tolist(),tiles=sorted(tiles),rowtiles=sorted(rowtiles)),ids[order],zs[order]
def main():
 common=np.load(ROOT/'galaxy_angular_support.npy')
 for i in range(2,7):common &=np.load(t.B/f'ph{i:03d}/p3_fields/angular_support_nside256.npz')['support']
 paths={'real_kibo':BASE/'LSS/kibo-v1/LSScats/v1/BGS_BRIGHT_full_HPmapcut.dat.fits','real_loa':BASE/'LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits','mock_kibo':BASE/'mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl6/kibo-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.fits','mock_loa':BASE/'mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl6/loa-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.h5'};results={};arr={}
 for key,path in paths.items():
  r,ids,z=scan(path,key.startswith('real'),common);results[key]=r;arr[key]=(ids,z);print(key,r['rows'],r['stages']['success'],len(r['tiles']),flush=True)
 joins={}
 for label,a,b in [('real','real_kibo','real_loa'),('mock','mock_kibo','mock_loa')]:
  ids,ia,ib=np.intersect1d(arr[a][0],arr[b][0],return_indices=True);dz=abs(arr[a][1][ia]-arr[b][1][ib]);joins[label]=dict(shared_successful_ids=len(ids),only_kibo=len(arr[a][0])-len(ids),only_loa=len(arr[b][0])-len(ids),redshift_changed_gt_0p001=int((dz>.001).sum()),tile_symmetric_difference=len(set(results[a]['tiles'])^set(results[b]['tiles'])))
 t.save(OUT/'PRODUCT_CROSSWALK.json',dict(products=results,joins=joins,job=os.environ.get('SLURM_JOB_ID'),script_sha256=t.digest(__file__),note='Real comparison includes LSS v1 versus v2.1 changes, not pure redshift pipeline difference. Same fixed common mask for counts; no catalogue changed.'));print(joins,flush=True)
if __name__=='__main__':main()
