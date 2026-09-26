"""Compare observed g-r at fixed fine z and numerical r, common southern sky.

Histogram standardization is diagnostic only; no catalogue/selection is changed.
"""
from pathlib import Path
import os,json,hashlib
import numpy as np,fitsio,healpy as hp
import loa_trial as t
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_photometry_recipe_20260926'
def main():
 root=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1');common=np.load(root/'galaxy_angular_support.npy')
 for i in range(2,7):common &=np.load(t.B/f'ph{i:03d}/p3_fields/angular_support_nside256.npz')['support']
 sources={'loa':t.read(root/'INPUTS_READY.json')['source']['path'],'ph006':str(t.B/'ph006/catalogues/observed/ph006_bgs_bright_full_observed_with_tweb.fits')}
 ze=np.linspace(.15,.55,41);re=np.linspace(12,19.5,76);ce=np.linspace(-2,5,351);hist={};excluded={}
 for name,path in sources.items():
  hist[name]=np.zeros((40,75,350),dtype=np.int64);excluded[name]=0
  cols=['RA','DEC','ZWARN']+(['Z_not4clus','DELTACHI2','SPECTYPE','PHOTSYS','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R'] if name=='loa' else ['Z','R_MAG_APP','G_R_OBS'])
  with fitsio.FITS(path) as f:
   for start in range(0,f[1].get_nrows(),250000):
    d=f[1].read(rows=np.arange(start,min(start+250000,f[1].get_nrows())),columns=cols);z=d['Z_not4clus'] if name=='loa' else d['Z'];k=(z>=.15)&(z<.55)&(d['ZWARN']==0)&(t.galactic_cap(d['RA'],d['DEC'])==0)&common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
    if name=='loa':
     k &=(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')&(np.char.strip(d['PHOTSYS'].astype('U'))=='S')
     r=22.5-2.5*np.log10(np.maximum(d['FLUX_R']/d['MW_TRANSMISSION_R'],1e-20));c=-2.5*np.log10(np.maximum(d['FLUX_G']/d['MW_TRANSMISSION_G'],1e-20))-(r-22.5)
    else:r=d['R_MAG_APP'];c=d['G_R_OBS']
    k &=(r>=12)&(r<19.5)&np.isfinite(c);excluded[name]+=int(np.count_nonzero(k&((c<ce[0])|(c>=ce[-1]))))
    hist[name]+=np.histogramdd(np.column_stack([z[k],r[k],c[k]]),bins=[ze,re,ce])[0].astype('i8')
  print(name,int(hist[name].sum()),flush=True)
 rows=[];cen=(ce[:-1]+ce[1:])/2
 for j in range(4):
  a=hist['loa'][j*10:(j+1)*10];b=hist['ph006'][j*10:(j+1)*10];na=a.sum(-1);nb=b.sum(-1);ok=(na>=20)&(nb>=20)
  w=np.divide(na,nb,out=np.zeros_like(na,dtype=float),where=ok)
  aa=a[ok].sum(0);bb=(b*w[...,None]).sum((0,1))
  med=lambda h:float(np.interp(.5*h.sum(),np.cumsum(h),cen))
  rows.append(dict(z=[float(ze[j*10]),float(ze[(j+1)*10])],cells=int(ok.sum()),loa_retained=int(aa.sum()),loa_total=int(na.sum()),mock_retained=int(b[ok].sum()),mock_total=int(nb.sum()),loa_median=med(aa),mock_standardized_median=med(bb),median_difference=med(aa)-med(bb),loa_mean=float(aa@cen/aa.sum()),mock_standardized_mean=float(bb@cen/bb.sum())))
 np.savez_compressed(OUT/'MATCHED_CELLS.npz',**hist,z_edges=ze,r_edges=re,colour_edges=ce)
 out=dict(rows=rows,sources=sources,excluded_colour_outside_hist=excluded,job=os.environ.get('SLURM_JOB_ID'),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),method='SGC only; Loa PHOTSYS=S; common angular support; both 12<=r<19.5 numerically. Cells dz=.01, dr=.1, >=20 galaxies per catalogue; mock cell colour distributions standardized to real cell counts. Passbands and magnitude estimators not assumed equivalent; no production changes.')
 (OUT/'MATCHED_CELLS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
