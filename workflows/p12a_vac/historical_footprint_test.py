"""Historical BRIGHT geometry test; diagnostic counts and geometric volumes only."""
import os,json,hashlib
from pathlib import Path
import numpy as np,healpy as hp
from astropy.table import Table
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from desimodel.footprint import is_point_in_desi
from alignment_parent_screen import BASE,MASK
from mock_selection_test import cap
R=Path(__file__).resolve().parents[2];O=R/'docs/evidence/p12a_historical_footprint_20260928';O.mkdir(exist_ok=True)
tp=Path('/global/common/software/desi/perlmutter/desiconda/20230111-2.1.0/code/desimodel/main/data/footprint/desi-tiles.ecsv')
t=Table.read(tp);dtype=str(t['IN_DESI'].dtype);t=t[(t['PROGRAM']=='BRIGHT')&(t['IN_DESI']==1)]
common=np.load(MASK);pix=np.flatnonzero(common);ra,dec=hp.pix2ang(256,pix,lonlat=True)
inside=np.zeros(len(common),bool);inside[pix]=is_point_in_desi(t,ra,dec)
interior=inside.copy()
for _ in range(2):
 ix=np.flatnonzero(interior);nn=hp.get_all_neighbours(256,ix);keep=np.all((nn>=0)&interior[np.maximum(nn,0)],axis=0);interior[ix[~keep]]=False
# Equal-area quadrature at two resolutions, preserve exact object geometry below.
area={}
for ns in [512,1024]:
 sums=np.zeros((3,2));step=500000
 for a in range(0,hp.nside2npix(ns),step):
  pp=np.arange(a,min(a+step,hp.nside2npix(ns)));rr,dd=hp.pix2ang(ns,pp,lonlat=True);p=hp.ang2pix(256,rr,dd,lonlat=True);ok=common[p];rr=rr[ok];dd=dd[ok];p=p[ok];c=cap(rr,dd)
  hist=is_point_in_desi(t,rr,dd)
  for j,s in enumerate([np.ones(len(p),bool),hist,hist&interior[p]]):
   for i in [0,1]:sums[j,i]+=np.count_nonzero(s&(c==i))*hp.nside2pixarea(ns,degrees=True)
 area[str(ns)]=sums.tolist()
 print('AREA',ns,sums,flush=True)
ze=np.linspace(.15,.55,41);paths={v:BASE/v/'z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits' for v in ['v0.1','v1']}
paths['loa']=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
counts={};meta={};stride=20
for name,path in paths.items():
 before=[path.stat().st_size,path.stat().st_mtime_ns];h=np.zeros((3,2,40));raw=np.zeros_like(h);flag=np.zeros((2,2),int)
 with fits.open(path,memmap=True) as f:
  cols=['RA','DEC']+(['Z_not4clus','ZWARN','DELTACHI2','SPECTYPE','FLUX_R','MW_TRANSMISSION_R','FRACZ_TILELOCID','FRAC_TLOBS_TILES','WEIGHT_ZFAIL'] if name=='loa' else ['Z','R_MAG_APP'])
  hasflag='IN_Y5' in f[1].columns.names
  if hasflag:cols+=['IN_Y5']
  for start in range(0,len(f[1].data),2000000):
   d={c:np.array(f[1].data[c][start:start+2000000:stride]) for c in cols};z=d['Z_not4clus'] if name=='loa' else d['Z'];p=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True);ok=common[p]&(z>=.15)&(z<.55)
   if name=='loa':
    ok&=(d['ZWARN']==0)&(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
    with np.errstate(invalid='ignore',divide='ignore'):r=22.5-2.5*np.log10(d['FLUX_R']/d['MW_TRANSMISSION_R']);w=d['WEIGHT_ZFAIL']/(d['FRACZ_TILELOCID']*d['FRAC_TLOBS_TILES'])
   else:r=d['R_MAG_APP'];w=np.ones(len(z))
   ok&=(r>=12)&(r<19.5);assert np.all(np.isfinite(w[ok])&(w[ok]>0))
   zz=z[ok];ww=w[ok];rr=d['RA'][ok];dd=d['DEC'][ok];p=p[ok];c=cap(rr,dd);historical=is_point_in_desi(t,rr,dd)
   if hasflag:
    ff=d['IN_Y5'][ok]==1
    for x in [0,1]:
     for y in [0,1]:flag[x,y]+=np.count_nonzero((ff==x)&(historical==y))
   for j,s in enumerate([np.ones(len(p),bool),historical,historical&interior[p]]):
    for i in [0,1]:
     take=s&(c==i);h[j,i]+=np.histogram(zz[take],ze,weights=ww[take]*stride)[0];raw[j,i]+=np.histogram(zz[take],ze)[0]
 assert before==[path.stat().st_size,path.stat().st_mtime_ns]
 counts[name]=h;meta[name]={'path':str(path),'stat':before,'IN_Y5_available':hasflag,'flag_vs_geometry':flag.tolist(),'sample_counts':raw.tolist()};print('DONE',name,flush=True)
grid=np.linspace(0,1,100001);om=(.02237+.1200+.00064420)/.6736**2;chi=cumulative_trapezoid(299792.458/100/np.sqrt(om*(1+grid)**3+1-om),grid,initial=0);vol=np.array(area['1024'])[:,:,None]*(np.pi/180)**2*np.diff(np.interp(ze,grid,chi)**3)/3
assert np.all(vol>0)
result={'area_deg2':area,'domains':['existing_common','common_and_historical_exact','common_and_historical_interior_two_pixel_rings'],'caps':['SGC','NGC'],'z_edges':ze.tolist(),'counts':{k:v.tolist() for k,v in counts.items()},'density_h3_Mpc3':{k:(v/vol).tolist() for k,v in counts.items()},'shell_ratios':{k:(v.reshape(3,2,4,10).sum(-1)/counts['loa'].reshape(3,2,4,10).sum(-1)).tolist() for k,v in counts.items() if k!='loa'},'sources':meta,'tiles_sha256':hashlib.sha256(tp.read_bytes()).hexdigest(),'tiles_path':str(tp),'tiles_IN_DESI_dtype':dtype,'tiles_selected':len(t),'mask_sha256':hashlib.sha256(MASK.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'job':os.environ.get('SLURM_JOB_ID'),'stride':stride,'limits':'Raw parents versus assignment*zfail corrected Loa; geometric common volumes, not official random-derived imaging/veto effective volumes. Pixel-centre quadrature; exact per-object historical tile membership. Interior erodes both common and historical pixel mask twice, then exact membership. Deterministic sample one exposed phase; no cosmic variance significance. No repairs.'}
(O/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['shell_ratios']))
