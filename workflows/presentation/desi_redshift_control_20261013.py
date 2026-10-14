#!/usr/bin/env python3
"""Matched-row raw/redshift/fine-mass+redshift check for the user's flatness question."""
from pathlib import Path
import argparse,json
import numpy as np
import fitsio,healpy as hp
import matplotlib.pyplot as plt
from PIL import Image
from desi_mass_movie_20261013 import PRODUCT,MASS_EDGES,MASS_COLORS,CLASS_NAMES,estimate,sha,dump,style

def main(out,render_only=False):
 if not render_only:
  product=json.loads(PRODUCT.read_text());assert sha(product['path'])==product['sha256']
  cols=['TARGETID','RA','DEC','Z','SUPPORTED','P_WEB','CIGALE_MATCHED','CIGALE_SPECTYPE','MASS_CG_15','SFR_CG_15']
  d=fitsio.read(product['path'],columns=cols);mass=d['MASS_CG_15'];sfr=d['SFR_CG_15'];z=d['Z'];k=d['SUPPORTED']&d['CIGALE_MATCHED']&(np.char.strip(d['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')&(z>=.2)&(z<.3)&np.isfinite(mass)&np.isfinite(sfr)&(mass>0)&(sfr>0)
  d=d[k];lm=np.log10(mass[k]);y=(np.log10(sfr[k])-lm<-11).astype(float);p=d['P_WEB'].astype(float);reg=hp.ang2pix(8,d['RA'],d['DEC'],lonlat=True);zc=np.minimum(3,np.floor((d['Z']-.2)/.025).astype(int));cells=np.floor((lm-8+1e-10)/.1).astype(int)*4+zc
  result=[]
  for lo,hi in zip(MASS_EDGES[:-1],MASS_EDGES[1:]):
   k=(lm>=lo)&(lm<hi);initial=estimate(y[k],p[k],cells[k],reg[k]);assert initial
   retained=k&np.isin(cells,initial['retained_cell_ids']);rr={}
   for name,c in [('raw',np.zeros(len(d),int)),('redshift',zc),('mass_redshift',cells)]:
    rr[name]=estimate(y[retained],p[retained],c[retained],reg[retained]);assert rr[name]['retained_rows']==int(retained.sum())
   result.append({'mass_range':[lo,hi],'identical_rows':int(retained.sum()),'estimators':rr})
  old=json.loads((out/'mass_statistics.json').read_text())
  for a,b in zip(result,old['curves']):assert np.allclose(a['estimators']['mass_redshift']['mean'],b['class']['mean'],atol=1e-12)
  dump(out/'redshift_control_statistics.json',{'curves':result,'source_sha256':product['sha256'],'script_sha256':sha(__file__),'reference':'Identical galaxies within each mass band for raw, common-redshift, and common-fine-mass/redshift estimators. Each controlled estimator uses its own fixed common cell mixture.','intervals':'16–84% of96 paired-within-estimator nside8 sky-block bootstrap draws. Contrasts computed within each replicate, not by subtracting marginal intervals. Systematics excluded.'})
 else:result=json.loads((out/'redshift_control_statistics.json').read_text())['curves']
 style();f,axs=plt.subplots(1,3,figsize=(16,9),dpi=120,sharey=True);f.subplots_adjust(left=.08,right=.98,bottom=.24,top=.67,wspace=.15)
 for j,(name,title) in enumerate([('raw','Raw within mass bands'),('redshift','Same redshift mix'),('mass_redshift','Same mass + redshift mix')]):
  ax=axs[j]
  for i,r in enumerate(result):
   v=r['estimators'][name];mu=np.array(v['mean']);ci=np.array(v['interval']);ax.plot(range(4),mu,color=MASS_COLORS[i],marker=['o','s','^','D'][i],ms=6,lw=2,label='–'.join(f'{n:.1f}' for n in r['mass_range']));ax.vlines(range(4),ci[0],ci[1],color=MASS_COLORS[i],lw=1)
  ax.set_xticks(range(4),CLASS_NAMES);ax.tick_params(axis='x',labelsize=12);ax.set_ylim(0,1);ax.set_title(title,fontsize=19,pad=17);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_xlabel('Inferred environment',fontsize=17,labelpad=12)
 axs[0].set_ylabel('Quenched fraction',fontsize=20);f.text(.5,.90,'Does redshift balancing change the trend?',fontsize=26,ha='center');f.text(.5,.84,'Identical galaxies in each panel · DESI + CIGALE · 0.20 ≤ z < 0.30',fontsize=17,ha='center');handles,labels=axs[0].get_legend_handles_labels();f.legend(handles,labels,ncol=4,frameon=False,loc='center',bbox_to_anchor=(.5,.78),fontsize=16,title='log₁₀(M★ / M☉)',title_fontsize=14)
 f.text(.5,.14,'Probability-weighted environments · quenched: log₁₀(sSFR / yr⁻¹) < −11',ha='center',fontsize=14,color='#C8CDD4');f.text(.5,.09,'Intervals: 16–84% sky-block bootstrap · real-DESI calibration unverified',ha='center',fontsize=14,color='#C8CDD4')
 for ext in ['svg','pdf','png']:f.savefig(out/f'redshift_control_check.{ext}',dpi=240 if ext=='png' else 120)
 plt.close(f);Image.open(out/'redshift_control_check.png').resize((1280,720),Image.Resampling.LANCZOS).save(out/'redshift_control_check_preview.png')
 print(json.dumps([{'mass':r['mass_range'],'contrasts_pp':{k:100*v['last_minus_first'] for k,v in r['estimators'].items()}} for r in result]),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);p.add_argument('--render-only',action='store_true');a=p.parse_args();main(a.out,a.render_only)
