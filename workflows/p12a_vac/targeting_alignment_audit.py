"""Audit observable BGS targeting components; do not claim full target replay."""
from pathlib import Path
import json,hashlib,os
import numpy as np, fitsio,healpy as hp
import loa_trial as t
ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_targeting_alignment_20260925.json'
def main():
 common=np.load(ROOT/'galaxy_angular_support.npy')
 for i in range(2,7):common &= np.load(t.B/f'ph{i:03d}/p3_fields/angular_support_nside256.npz')['support']
 sources={'loa':t.read(ROOT/'INPUTS_READY.json')['source']['path'],'ph006':str(t.B/'ph006/catalogues/observed/ph006_bgs_bright_full_observed_with_tweb.fits')};result={};edges=np.linspace(-2,5,141);medges=np.linspace(12,21,181)
 for kind,path in sources.items():
  counts={};grhist=np.zeros((8,140),dtype='i8');mhist=np.zeros((8,180),dtype='i8');headers=[];finehist=np.zeros((2,40,140),dtype='i8')
  with fitsio.FITS(path) as f:
   columns=f[1].get_colnames()
   for hdu in f:
    headers.append([{k:r[k] for k in ['name','value','comment'] if k in r} for r in hdu.read_header().records() if not r['name'].startswith(('TTYPE','TFORM','TUNIT','TDIM'))])
   for start in range(0,f[1].get_nrows(),250000):
    d=f[1][start:min(start+250000,f[1].get_nrows())];z=d['Z_not4clus'] if kind=='loa' else d['Z'];shell=np.searchsorted([.15,.25,.35,.45,.55],z,side='right')-1;k=np.isfinite(z)&(shell>=0)&(shell<4)&(d['ZWARN']==0)
    if kind=='loa':k&=(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
    d=d[k];shell=shell[k];pix=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True);k=common[pix];d=d[k];shell=shell[k];cap=t.galactic_cap(d['RA'],d['DEC']);b=cap*4+shell
    checks={'baseline':np.ones(len(d),bool),'bright_bit':(d['BGS_TARGET']&2)!=0}
    if kind=='loa':
     flux={c:d['FLUX_'+c]/d['MW_TRANSMISSION_'+c] for c in ['G','R','Z']};mag=lambda v:22.5-2.5*np.log10(np.maximum(v,1e-16));r=mag(flux['R']);gr=mag(flux['G'])-r;rfib=mag(d['FIBERFLUX_R']/d['MW_TRANSMISSION_R']);north=np.char.strip(d['PHOTSYS'].astype('U'))=='N';south=np.char.strip(d['PHOTSYS'].astype('U'))=='S'
     checks.update(photsys_N=north,photsys_S=south,nominal_bright_magnitude=(r>=12)&(r<np.where(north,19.54,19.5)),uniform_r19p5=r<19.5,fiber_magnitude=((rfib<5.1+r)&(r<17.8))|((rfib<22.9)&(r>17.8)&(r<20))|((rfib<2.9+r)&(r>20)),fiber_total_limit=d['FIBERTOTFLUX_R']<=10**((22.5-15)/2.5),positive_ivar=(d['FLUX_IVAR_G']>0)&(d['FLUX_IVAR_R']>0)&(d['FLUX_IVAR_Z']>0),colour_box=(flux['R']>flux['G']*10**(-1/2.5))&(flux['R']<flux['G']*10**(4/2.5))&(flux['Z']>flux['R']*10**(-1/2.5))&(flux['Z']<flux['R']*10**(4/2.5)),imaging_bits_1_13=(d['MASKBITS']&((1<<1)|(1<<13)))==0,maskbit11_clear=(d['MASKBITS']&(1<<11))==0,nobs=(d['NOBS_G']>0)&(d['NOBS_R']>0)&(d['NOBS_Z']>0),dchi40=d['DELTACHI2']>40,tsnr1000=d['TSNR2_BGS']>1000)
     checks['available_target_components']=np.logical_and.reduce([checks[n] for n in ['nominal_bright_magnitude','fiber_magnitude','fiber_total_limit','positive_ivar','colour_box','imaging_bits_1_13','nobs']])
     checks['strict_joint_diagnostic']=checks['available_target_components']&checks['uniform_r19p5']&checks['dchi40']
    else:r=d['R_MAG_APP'];gr=d['G_R_OBS'];checks['uniform_r19p5']=r<19.5
    for label,k in checks.items():counts.setdefault(label,np.zeros(8,dtype='i8'));counts[label]+=np.bincount(b[k],minlength=8)
    fine=np.searchsorted(np.linspace(.15,.55,41),d['Z_not4clus'] if kind=='loa' else d['Z'],side='right')-1
    for cc in [0,1]:
     for jj in range(40):finehist[cc,jj]+=np.histogram(gr[(cap==cc)&(fine==jj)&checks['uniform_r19p5']],edges)[0]
    for j in range(8):grhist[j]+=np.histogram(gr[b==j],edges)[0];mhist[j]+=np.histogram(r[b==j],medges)[0]
  quant=lambda h,e:[[float(np.interp(q*np.sum(row),np.cumsum(row),(e[1:]+e[:-1])/2)) for q in [.16,.5,.84]] for row in h]
  result[kind]=dict(path=path,columns=columns,headers=headers,counts={k:v.tolist() for k,v in counts.items()},fine_gr_hist_uniform_r19p5=finehist.tolist(),gr_hist=grhist.tolist(),r_hist=mhist.tolist(),gr_quantiles_hist=quant(grhist,edges),r_quantiles_hist=quant(mhist,medges))
 paths=[Path(__file__),Path('/global/common/software/desi/perlmutter/desiconda/current/code/desitarget/main/py/desitarget/cuts.py'),Path('/global/u2/d/dkololgi/TNG/Illustris/workflows/abacus_tweb/secondgen_mocks/ph000/scripts/upstream_prepare_mocks_Y3_bright.py')]
 out=dict(job=os.environ.get('SLURM_JOB_ID'),catalogues=result,gr_edges=edges.tolist(),r_edges=medges.tolist(),sources_sha256={str(p):t.digest(p) for p in paths},order='SGC shells0–3 then NGC shells0–3',full_selection_alignment_verified=False,note='Current installed main targeting formulas are diagnostic, not proven historical production version. Loa table lacks REF_CAT and Gaia inputs: no full isBGS replay or SGA recovery. FIBERTOTFLUX is uncorrected as in desitarget. Individual tests need not all hold for legitimate SGA recovery. No changes to selection or inference.')
 t.save(OUT,out);print(json.dumps({k:v['counts'] for k,v in result.items()},indent=2),flush=True)
if __name__=='__main__':main()
