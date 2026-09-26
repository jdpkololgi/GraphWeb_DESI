"""Full-survey descriptive posterior diagnostics; no DESI truth/calibration claims."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import fitsio,healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def main(root,out):
 out.mkdir(parents=True,exist_ok=False)
 complete=json.loads((root/'FULL_VAC_COMPLETE.json').read_text());path=Path(complete['catalogue']);assert digest(path)==complete['catalogue_sha256']
 d=fitsio.read(path);yes=d['SUPPORTED'];v=d[yes];del d
 assert len(v)==complete['supported_rows'];assert np.isfinite(v['P_WEB']).all();assert np.allclose(v['P_WEB'].sum(axis=1),1)
 med=v['EIGENVALUE_Q50'];width=v['EIGENVALUE_Q84']-v['EIGENVALUE_Q16'];assert np.all(width>=0)
 z=v['Z'];p=v['P_WEB'];conf=p.max(axis=1);colors=['#2166ac','#c07818','#aa3377'];names=['Void','Sheet','Filament','Knot'];cc=['#4477aa','#ccbb44','#228833','#aa3377']
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#252525','text.color':'#252525','figure.facecolor':'white','savefig.facecolor':'white','axes.grid':True,'grid.alpha':.16})
 pdf=PdfPages(out/'DESI_posterior_diagnostics.pdf');files=[]
 def save(fig,name):
  fig.savefig(out/(name+'.png'),dpi=170,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');files.append(name+'.png');plt.close(fig)
 bins=np.linspace(.15,.55,17);zc=(bins[:-1]+bins[1:])/2;wq=[];pm=[];counts=[]
 for a,b in zip(bins[:-1],bins[1:]):
  k=(z>=a)&(z<b);counts.append(int(k.sum()));wq.append(np.quantile(width[k],[.16,.5,.84],axis=0));pm.append(p[k].mean(axis=0,dtype='f8'))
 wq=np.array(wq);pm=np.array(pm)
 fig,ax=plt.subplots(2,3,figsize=(14,8),layout='constrained');fig.suptitle('DESI Loa environmental posteriors — full survey, provisional',fontsize=18)
 for j in range(3):
  a=ax[0,j];bounds=np.quantile(med[:,j],[.001,.999]);edges=np.linspace(*bounds,100)
  for lo,hi,style,col in [(.15,.35,'-','#2166ac'),(.35,.55,'--','#b36716')]:
   x=med[(z>=lo)&(z<hi),j];h,_=np.histogram(x,edges);a.stairs(h/(len(x)*np.diff(edges)),edges,label=f'{lo:.2f} ≤ z < {hi:.2f}',color=col,linestyle=style,lw=1.6)
  a.axvline(.2,color='#555555',ls=':',lw=1);a.set(xlabel=fr'Posterior median $\lambda_{j+1}$',ylabel='Galaxy density per eigenvalue unit',title=fr'$\lambda_{j+1}$: distribution of galaxy medians');a.legend(fontsize=9)
 a=ax[1,0]
 for j in range(3):a.plot(zc,wq[:,1,j],color=colors[j],marker=['o','s','^'][j],label=fr'$\lambda_{j+1}$');a.fill_between(zc,wq[:,0,j],wq[:,2,j],color=colors[j],alpha=.10)
 a.axvspan(.35,.55,color='#b36716',alpha=.08);a.set(xlabel='Observed redshift',ylabel='Q84 − Q16',title='Width of each galaxy’s central 68% interval');a.legend();a.text(.02,.96,'Bands: 16–84% spread across galaxies',transform=a.transAxes,va='top',fontsize=9)
 a=ax[1,1]
 for j in range(4):a.plot(zc,pm[:,j],label=names[j],color=cc[j],marker=['o','s','^','d'][j])
 a.axvspan(.35,.55,color='#b36716',alpha=.08);a.set(xlabel='Observed redshift',ylabel='Mean class probability per galaxy',ylim=(0,1),title='Posterior environment probabilities');a.legend(fontsize=9,ncol=2)
 a=ax[1,2]
 for lo,hi,style,col in [(.15,.35,'-','#2166ac'),(.35,.55,'--','#b36716')]:
  k=(z>=lo)&(z<hi);a.hist(conf[k],bins=np.linspace(.25,1,61),weights=np.full(k.sum(),1/k.sum()),histtype='step',color=col,ls=style,lw=1.6,label=f'{lo:.2f} ≤ z < {hi:.2f}')
 a.set(xlabel='Largest of the four class probabilities',ylabel='Fraction of galaxies per bin',title='How concentrated are the class probabilities?');a.legend(fontsize=9)
 fig.supxlabel(f'{len(v):,} supported galaxies · 512 draws each · threshold λ = 0.2 · R = 7 Mpc/h, target epoch z = 0.2\nShaded/high-z group: unresolved selection mismatch. Descriptive diagnostics, not real-data coverage.',fontsize=10)
 save(fig,'01_posterior_overview')
 # Equal-area HEALPix sky pixels, deliberately restrict redshift to reduce known high-z mismatch.
 k=(z>=.15)&(z<.25);pix=hp.ang2pix(64,v['RA'][k],v['DEC'][k],lonlat=True);n=np.bincount(pix,minlength=hp.nside2npix(64));valid=n>=10;pp=np.flatnonzero(valid);ra,dec=hp.pix2ang(64,pp,lonlat=True)
 fig=plt.figure(figsize=(14,5.7),layout='constrained');fig.suptitle('Sky projection — galaxies at 0.15 ≤ z < 0.25',fontsize=17)
 for j,(values,title,cmap) in enumerate([(p[k,2],'Mean filament probability','viridis'),(width[k,1],'Mean 68% interval width of λ₂','magma')]):
  mean=np.bincount(pix,weights=values,minlength=len(n))[valid]/n[valid];a=fig.add_subplot(1,2,j+1,projection='mollweide');im=a.scatter(-np.deg2rad(ra-180),np.deg2rad(dec),c=mean,s=8,cmap=cmap,rasterized=True,vmin=0 if j==0 else None,vmax=1 if j==0 else None);a.set_xticks(np.deg2rad([-120,-60,0,60,120]),labels=['300°','240°','180°','120°','60°']);a.set_title(title);a.set_xlabel('ICRS right ascension');a.set_ylabel('Declination');fig.colorbar(im,ax=a,orientation='horizontal',shrink=.8,pad=.12)
 fig.supxlabel('Per-galaxy averages in nside=64 equal-area pixels with ≥10 supported galaxies. Radial projection, not a reconstructed 3D field.\nProvisional model-conditional estimates; angular patterns can reflect selection and response.',fontsize=10);save(fig,'02_sky_projection')
 # Four reproducibly selected individual examples, one per shell, typical lambda2 width away from boundaries.
 chosen=[]
 for lo,hi in zip([.15,.25,.35,.45],[.25,.35,.45,.55]):
  ii=np.flatnonzero((z>=lo)&(z<hi)&(v['BOUNDARY_MPC']>20.6917)&((v['QUALITY']&16)==0));target=np.median(width[ii,1]);chosen.append(int(ii[np.argmin(np.abs(width[ii,1]-target))]))
 index={a['core_id']:a for a in json.loads((root/'DRAW_INDEX.json').read_text())['shards']};draws=[];examples=[]
 for row in chosen:
  item=index[int(v['CORE_ID'][row])];assert digest(item['path'])==item['sha256']
  with np.load(item['path']) as f:
   ix=np.flatnonzero(f['TARGETID']==v['TARGETID'][row]);assert len(ix)==1;draws.append(f['eigenvalue_draws'][ix[0]])
  examples.append(dict(targetid=int(v['TARGETID'][row]),redshift=float(z[row]),quality=int(v['QUALITY'][row]),p_web=p[row].tolist(),source=item['path']))
 fig,ax=plt.subplots(4,3,figsize=(13,11),layout='constrained');fig.suptitle('Individual galaxies — saved posterior draws',fontsize=18)
 for i,sample in enumerate(draws):
  for j in range(3):
   a=ax[i,j];a.hist(sample[:,j],bins=32,density=True,histtype='stepfilled',color=colors[j],alpha=.35);a.axvline(.2,color='#555555',ls=':',lw=1);a.axvline(med[chosen[i],j],color=colors[j],lw=1.4);a.set_xlabel(fr'$\lambda_{j+1}$');a.set_ylabel('Posterior density')
   if j==1:a.set_title(f"z={examples[i]['redshift']:.3f} · TARGETID {examples[i]['targetid']}",fontsize=10)
   if j==2:a.text(.97,.95,'P(void, sheet, filament, knot)\n'+', '.join(f'{q:.2f}' for q in examples[i]['p_web']),transform=a.transAxes,ha='right',va='top',fontsize=8)
 fig.supxlabel('One galaxy per redshift shell, chosen nearest the median λ₂ width among supported interior galaxies. Not a random sample.\n512 saved draws each; vertical dotted line is class threshold 0.2. Last two shells have unresolved selection mismatch.',fontsize=10);save(fig,'03_individual_posteriors')
 fig,ax=plt.subplots(2,3,figsize=(12,7),layout='constrained');fig.suptitle('Joint eigenvalue posteriors — two of the selected examples',fontsize=17)
 for i,which in enumerate([0,3]):
  sample=draws[which]
  for j,(a,b) in enumerate([(0,1),(0,2),(1,2)]):
   q=ax[i,j];q.scatter(sample[:,a],sample[:,b],s=5,alpha=.3,color=colors[i]);q.axvline(.2,color='grey',ls=':');q.axhline(.2,color='grey',ls=':');q.set(xlabel=fr'$\lambda_{a+1}$',ylabel=fr'$\lambda_{b+1}$',title=f"z={examples[which]['redshift']:.3f}")
 fig.supxlabel('Every point is one joint draw; marginal posteriors are dependent. The high-z example is selection-flagged.\nNo DESI ground-truth eigenvalues are available for a coverage test.',fontsize=10);save(fig,'04_joint_posteriors');pdf.close()
 summary=dict(rows=complete['rows'],supported=len(v),unsupported=complete['rows']-len(v),catalogue_sha256=complete['catalogue_sha256'],mean_class_probabilities=p.mean(axis=0,dtype='f8').tolist(),median_max_class_probability=float(np.median(conf)),fraction_max_probability_ge_0p8=float(np.mean(conf>=.8)),median_68_width=np.median(width,axis=0).tolist(),redshift_centres=zc.tolist(),bin_counts=counts,width_population_quantiles=wq.tolist(),class_probabilities_by_redshift=pm.tolist(),examples=examples,figures=files,script_sha256=digest(__file__),science_release_ready=False)
 (out/'PLOT_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:summary[k] for k in ['supported','mean_class_probabilities','median_max_class_probability','fraction_max_probability_ge_0p8','median_68_width']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);a=p.parse_args();main(a.root,a.output)
