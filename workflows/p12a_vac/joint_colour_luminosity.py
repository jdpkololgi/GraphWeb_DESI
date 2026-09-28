"""Exploratory common-photometry comparison; deterministic stride, no catalogue changes."""
import sys,json,hashlib,os
from pathlib import Path
import numpy as np, fitsio, healpy as hp
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from alignment_parent_screen import BASE,MASK
R=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(R/'docs/evidence/p12a_canonical_v1_replay_20260928/source'))
from hodpy.k_correction import DESI_KCorrection,DESI_KCorrection_color
OUT=R/'docs/evidence/p12a_joint_colour_luminosity_20260928';OUT.mkdir(exist_ok=True)
ZE=np.linspace(.15,.55,41);RE=np.linspace(12,19.5,76);CE=np.linspace(-2,5,141);ME=np.linspace(-27,-14,131)
STRIDE=20
paths={v:BASE/v/'z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits' for v in ['v0.1','v1']}
paths['loa']=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
sky=np.load(MASK);grid=np.linspace(0,1,100001);om=(.02237+.1200+.00064420)/.6736**2
chi=cumulative_trapezoid(299792.458/100/np.sqrt(om*(1+grid)**3+1-om),grid,initial=0)
k={p:DESI_KCorrection('r',p) for p in ['N','S']};kc={p:DESI_KCorrection_color(p) for p in ['N','S']}
allh={};meta={}
for name,path in paths.items():
 h=np.zeros((2,40,75,140));hm=np.zeros((2,40,130));nsel=np.zeros(2,dtype=int);nout=np.zeros(2,dtype=int);mout=np.zeros(2,dtype=int)
 cols=['RA','DEC']+(['Z_not4clus','ZWARN','DELTACHI2','SPECTYPE','PHOTSYS','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R','FRACZ_TILELOCID','FRAC_TLOBS_TILES','WEIGHT_ZFAIL'] if name=='loa' else ['Z','R_MAG_APP','G_R_OBS'])
 stat=[path.stat().st_size,path.stat().st_mtime_ns]
 with fits.open(path,memmap=True) as f:
  for start in range(0,f[1].header['NAXIS2'],2000000):
   d={col:np.array(f[1].data[col][start:min(start+2000000,f[1].header['NAXIS2']):STRIDE]) for col in cols}
   z=d['Z_not4clus'] if name=='loa' else d['Z'];sel=(z>=.15)&(z<.55)&sky[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
   if name=='loa':
    sel &=(d['ZWARN']==0)&(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
    fg=d['FLUX_G']/d['MW_TRANSMISSION_G'];fr=d['FLUX_R']/d['MW_TRANSMISSION_R'];sel &=(fg>0)&(fr>0)
    with np.errstate(invalid='ignore',divide='ignore'):r=22.5-2.5*np.log10(fr);c=-2.5*np.log10(fg/fr)
    ps=np.char.strip(d['PHOTSYS'].astype('U'));w=d['WEIGHT_ZFAIL']/(d['FRACZ_TILELOCID']*d['FRAC_TLOBS_TILES'])
   else:r=d['R_MAG_APP'];c=d['G_R_OBS'];ps=np.where(d['DEC']>32.375,'N','S');w=np.ones(len(d['RA']))
   sel &=(r>=12)&(r<19.5)&np.isfinite(c)
   assert np.all(np.isfinite(w[sel])&(w[sel]>0))
   for i,p in enumerate(['S','N']):
    take=sel&(ps==p);zz=z[take];rr=r[take];cc=c[take];ww=w[take]
    nsel[i]+=len(zz);nout[i]+=np.count_nonzero((cc<CE[0])|(cc>=CE[-1]))
    h[i]+=np.histogramdd(np.column_stack([zz,rr,cc]),[ZE,RE,CE],weights=ww)[0]
    rest=kc[p].rest_frame_colour(zz,cc);dm=5*np.log10((1+zz)*np.interp(zz,grid,chi))+25
    mag=rr-dm-k[p].k(zz,rest)+.67*(zz-.1)
    mout[i]+=np.count_nonzero(~np.isfinite(mag)|(mag<ME[0])|(mag>=ME[-1]))
    hm[i]+=np.histogram2d(zz,mag,bins=[ZE,ME],weights=ww)[0]
   print(name,start,flush=True)
 assert stat==[path.stat().st_size,path.stat().st_mtime_ns]
 allh[name]=h;allh[name+'_M']=hm;meta[name]={'path':str(path),'stat':stat,'sample_selected':nsel.tolist(),'colour_outside':nout.tolist(),'magnitude_outside':mout.tolist()}
rows=[];cen=(CE[:-1]+CE[1:])/2
for v in ['v0.1','v1']:
 for i,p in enumerate(['S','N']):
  for j in range(4):
   a=allh['loa'][i,10*j:10*j+10];b=allh[v][i,10*j:10*j+10];na=a.sum(-1);nb=b.sum(-1);ok=(na>=20)&(nb>=20)
   aa=a[ok].sum(0);bb=(b*np.divide(na,nb,out=np.zeros_like(na),where=ok)[...,None]).sum((0,1))
   norm_a=aa/aa.sum();norm_b=bb/bb.sum()
   rows.append({'mock':v,'photsys':p,'zlo':float(ZE[10*j]),'zhi':float(ZE[10*j+10]),'count_ratio':float(nb.sum()/na.sum()),'matched_cells':int(ok.sum()),'loa_weight_fraction_retained':float(aa.sum()/na.sum()),'mean_colour_mock_minus_loa':float((norm_b-norm_a)@cen),'colour_total_variation':float(.5*np.abs(norm_b-norm_a).sum())})
np.savez_compressed(OUT/'HISTOGRAMS.npz',**allh,z_edges=ZE,r_edges=RE,colour_edges=CE,M_edges=ME)
(OUT/'RESULTS.json').write_text(json.dumps({'rows':rows,'sources':meta,'stride':STRIDE,'job':os.environ.get('SLURM_JOB_ID'),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'contract':'Common NSIDE256 support and numeric 12<=r<19.5; Loa assignment*zfail weights; PHOTSYS S/N, mocks DEC32.375. Raw parents vs corrected selected data. Conditional colour dz=.01 dr=.1 cells >=20 weighted counts, standardization is diagnostic only. M is a common v1 K/E mapping hypothesis with approximate native distance, not measured independent luminosity or a volume-corrected LF. Deterministic stride20 is exploratory, not independent statistical sampling.'},indent=2)+'\n')
print(json.dumps(rows,indent=2))
