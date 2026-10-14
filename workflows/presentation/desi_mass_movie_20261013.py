#!/usr/bin/env python3
"""Current-P12-A mass-controlled slide and actual-galaxy rotation; no inference.
Run on a compute node: python THIS --out SCRATCH_DIRECTORY [--stage all|prepare|render]
"""
from pathlib import Path
import argparse, hashlib, json, os, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.animation import FFMpegWriter
from PIL import Image
ROOT=Path(os.environ.get('GRAPHWEB_REPO_ROOT','/global/homes/d/dkololgi/GraphWeb_DESI')).resolve()
PRODUCT=ROOT/'docs/figures/p12a_loa_cigale_20260929/PRODUCT.json'
MASS_EDGES=np.array([9.8,10.4,10.6,11.,11.6])
TRACE_EDGES=np.array([-np.inf,-.1,.1,.3,.5,.7,1.,1.4,np.inf])
CLASS_NAMES=['Void','Wall','Filament','Cluster']
CLASS_COLORS=['#9CF6D5','#5B8FF9','#E85D9E','#F2C14E']
MASS_COLORS=['#BD69EA','#E85D9E','#FF984B','#F3D842']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def estimate(y,w,cells,regions,x=None,nboot=96):
 """Common fine-mass/redshift mix across G soft or disjoint groups."""
 cell_ids,cells=np.unique(cells,return_inverse=True);_,regions=np.unique(regions,return_inverse=True)
 nc=cells.max()+1;nr=regions.max()+1;ng=w.shape[1];flat=cells*nr+regions
 c=np.zeros((ng,nc,nr));s=c.copy();q=c.copy();sx=c.copy()
 for j in range(ng):
  for dest,weights in [(c,w[:,j]),(s,w[:,j]*y),(q,w[:,j]**2)]:
   dest[j]=np.bincount(flat,weights=weights,minlength=nc*nr).reshape(nc,nr)
  if x is not None:sx[j]=np.bincount(flat,weights=w[:,j]*x,minlength=nc*nr).reshape(nc,nr)
 n=c.sum(2);eff=n*n/np.maximum(q.sum(2),1e-30);common=(eff>=20).all(0)
 if not common.any():return None
 ref=n[:,common].sum(0);ref/=ref.sum()
 def calc(a,b):return (b[:,common]/a[:,common]*ref).sum(1)
 mean=calc(n,s.sum(2));rng=np.random.default_rng(20261013);boots=[]
 for _ in range(nboot):
  mult=np.bincount(rng.integers(nr,size=nr),minlength=nr);a=c@mult;b=s@mult
  if (a[:,common]>0).all():boots.append(calc(a,b))
 assert len(boots)>=.8*nboot
 out={'mean':mean.tolist(),'last_minus_first':float(mean[-1]-mean[0]),'last_minus_first_interval':np.quantile(np.asarray(boots)[:,-1]-np.asarray(boots)[:,0],[.16,.84]).tolist(),'interval':np.quantile(boots,[.16,.84],axis=0).tolist(),'input_rows':len(y),'retained_rows':int(common[cells].sum()),'cells':int(common.sum()),'sky_blocks':int(nr),'bootstrap_replicates':len(boots),'reference_weights':ref.tolist(),'retained_cell_ids':cell_ids[common].tolist()}
 if x is not None:out['x']=calc(n,sx.sum(2)).tolist()
 return out

def prepare(out):
 import fitsio,healpy as hp
 from astropy.cosmology import Planck18
 product=json.loads(PRODUCT.read_text());src=Path(product['path']);assert sha(src)==product['sha256']
 cols=['TARGETID','RA','DEC','Z','SUPPORTED','P_WEB','EIGENVALUE_MEAN','CIGALE_MATCHED','CIGALE_SPECTYPE','MASS_CG_15','SFR_CG_15']
 d=fitsio.read(src,columns=cols);assert len(d)==product['rows'];p=d['P_WEB'].astype(float);z=d['Z']
 supported=d['SUPPORTED'];assert np.isfinite(p[supported]).all() and np.allclose(p[supported].sum(1),1,atol=1e-6)
 mass=d['MASS_CG_15'];sfr=d['SFR_CG_15'];valid=supported&d['CIGALE_MATCHED']&(np.char.strip(d['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')&np.isfinite(mass)&np.isfinite(sfr)&(mass>0)&(sfr>0)&(z>=.2)&(z<.3)
 v=d[valid];lm=np.log10(mass[valid]);ss=np.log10(sfr[valid])-lm;tr=v['EIGENVALUE_MEAN'].astype(float).sum(1);assert np.isfinite(tr).all()
 reg=hp.ang2pix(8,v['RA'],v['DEC'],lonlat=True);cells=np.floor((lm-8+1e-10)/.1).astype(int)*4+np.minimum(3,np.floor((v['Z']-.2)/.025).astype(int))
 group=np.digitize(tr,TRACE_EDGES[1:-1]);curves=[]
 for lo,hi in zip(MASS_EDGES[:-1],MASS_EDGES[1:]):
  k=(lm>=lo)&(lm<hi);y=(ss[k]<-11).astype(float)
  a=estimate(y,v['P_WEB'][k].astype(float),cells[k],reg[k]);b=estimate(y,np.eye(8)[group[k]],cells[k],reg[k],x=tr[k])
  assert a is not None and b is not None,(lo,hi)
  curves.append({'mass_range':[lo,hi],'class':a,'trace':b})
 dump(out/'mass_statistics.json',{'mass_edges':MASS_EDGES.tolist(),'trace_edges':['-inf',-.1,.1,.3,.5,.7,1.,1.4,'inf'],'valid_property_rows':len(v),'curves':curves})
 # Geometrical wedge selected before reading probabilities; no CIGALE requirement.
 k=supported&(z>=.2)&(z<.3)&(d['RA']>=120)&(d['RA']<160)&(d['DEC']>=14.5)&(d['DEC']<30.6)
 rows=np.flatnonzero(k);eligible=len(rows);rng=np.random.default_rng(20261013)
 if len(rows)>18000:rows=np.sort(rng.choice(rows,18000,replace=False))
 w=d[rows];zz=np.linspace(.15,.55,10001);rr=np.interp(w['Z'],zz,Planck18.comoving_distance(zz).value)
 ra=np.deg2rad(w['RA']);dec=np.deg2rad(w['DEC']);xyz=rr[:,None]*np.column_stack([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
 a,b=np.deg2rad([140,22.55]);basis=np.array([[-np.sin(a),np.cos(a),0],[-np.sin(b)*np.cos(a),-np.sin(b)*np.sin(a),np.cos(b)],[np.cos(b)*np.cos(a),np.cos(b)*np.sin(a),np.sin(b)]])
 xyz=xyz@basis.T
 np.savez_compressed(out/'movie_galaxies.npz',TARGETID=w['TARGETID'],catalogue_rows=rows,RA=w['RA'],DEC=w['DEC'],Z=w['Z'],xyz_Mpc=xyz,P_WEB=w['P_WEB'])
 dump(out/'provenance.json',{'source':product,'product_receipt_sha256':sha(PRODUCT),'script_sha256':sha(__file__),'selection':'Supported unique CIGALE galaxy match, positive finite 15-band mass/SFR, 0.20 <= z < 0.30; full footprint. Four mass bins chosen to match supplied reference.','quenched_definition':'log10(sSFR/yr^-1) < -11','controls':'Each mass band independently: common 0.1-dex mass / 0.025-redshift cell mixture across all four P_WEB groups (left) or eight bins of posterior-mean trace (right). Minimum Kish effective count 20 per group per cell. Panels can retain different common support.','uncertainty':'16-84% percentile intervals of 96 nside=8 sky-block bootstrap replicates, seed20261013, fixed common reference mix. Not property-error, model-transfer or calibrated posterior uncertainty.','class_order':CLASS_NAMES,'trace_definition':'sum(EIGENVALUE_MEAN), expectation of ordered tidal eigenvalue sum; dimensionless. Binning the posterior mean, not a density truth measurement.','movie':{'selection':'Supported, RA[120,160), DEC[14.5,30.6), z[.20,.30); no CIGALE/property/probability cut.','eligible_rows':eligible,'display_rows':len(rows),'sample':'Uniform without replacement, seed20261013, sorted catalogue row order; independent of probabilities.','coordinate':'Planck18 observed-redshift comoving Cartesian coordinates, rotated into east/north/line-of-sight at RA140 DEC22.55, Mpc not Mpc/h. Distance uses interpolation on 10001-point z=.15-.55 grid.','color':'argmax(P_WEB) using Void/Wall/Filament/Cluster; alpha=.20+.80*max(P_WEB). White-to-color intro is presentation only. No posterior resampling or galaxy motion.','interpretation':'Rotation of fixed observed redshift-space galaxy positions with provisional per-galaxy inferred environments. Not cosmic time evolution or a spatially joint posterior field.'},'scientific_limits':'Provisional P12-A Loa, real-DESI calibration unverified. Associations are descriptive and not causal. z<.3 chosen to reduce known high-z population mismatch; does not establish validity.','no_new_inference':True,'slurm_job':os.getenv('SLURM_JOB_ID')})
 print('Prepared statistics and',len(rows),'movie galaxies',flush=True)

def style():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16,'text.color':'#F5F5F5','axes.labelcolor':'#F5F5F5','xtick.color':'#C8CDD4','ytick.color':'#C8CDD4','axes.edgecolor':'#899099','axes.facecolor':'black','figure.facecolor':'black','savefig.facecolor':'black','pdf.fonttype':42,'svg.fonttype':'none'})
def figure(out):
 style();r=json.loads((out/'mass_statistics.json').read_text());f,axs=plt.subplots(1,2,figsize=(16,9),dpi=120,sharey=True)
 f.subplots_adjust(left=.085,right=.97,bottom=.23,top=.69,wspace=.22)
 f.text(.5,.87,'Quenched fraction within stellar-mass bands',ha='center',fontsize=25)
 f.text(.5,.815,'DESI + CIGALE · 0.20 ≤ z < 0.30 · same mass and redshift mix within each band',ha='center',fontsize=16)
 for j,c in enumerate(r['curves']):
  lo,hi=c['mass_range'];color=MASS_COLORS[j]
  for ax,key in zip(axs,['class','trace']):
   v=c[key];x=np.arange(4) if key=='class' else np.array(v['x']);mu=np.array(v['mean']);ci=np.array(v['interval'])
   ax.plot(x,mu,color=color,marker=['o','s','^','D'][j],ms=7,lw=2.2,label=f'{lo:.1f}–{hi:.1f}')
   ax.vlines(x,ci[0],ci[1],color=color,lw=1.7);ax.fill_between(x,ci[0],ci[1],color=color,alpha=.13)
 for ax in axs:
  ax.set_ylim(0,1);ax.set_yticks([0,.25,.5,.75,1]);ax.grid(axis='y',alpha=.16);ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=16)
 axs[0].set_xticks(range(4),CLASS_NAMES);axs[0].set_xlabel('Inferred environment',fontsize=19,labelpad=12);axs[0].set_ylabel('Quenched fraction',fontsize=19,labelpad=12)
 axs[1].set_xlabel('Posterior mean tidal trace  E[λ₁ + λ₂ + λ₃]',fontsize=18,labelpad=12)
 axs[0].set_title('Probability-weighted classes',fontsize=19,pad=18);axs[1].set_title('Continuous environment',fontsize=19,pad=18)
 handles,labels=axs[0].get_legend_handles_labels();f.legend(handles,labels,ncol=4,loc='center',bbox_to_anchor=(.55,.753),frameon=False,fontsize=17,columnspacing=1.4)
 f.text(.095,.753,'log₁₀(M★/M☉)',fontsize=17,va='center')
 f.text(.5,.115,'Quenched: log₁₀(sSFR / yr⁻¹) < −11 · intervals: 16–84% sky-block bootstrap',ha='center',fontsize=14,color='#C8CDD4')
 f.text(.5,.075,'Provisional environments · descriptive associations; real-DESI calibration remains unverified',ha='center',fontsize=14,color='#C8CDD4')
 for ext in ['svg','pdf','png']:f.savefig(out/f'mass_control.{ext}',dpi=240 if ext=='png' else 120)
 plt.close(f);Image.open(out/'mass_control.png').resize((1280,720),Image.Resampling.LANCZOS).save(out/'mass_control_preview.png')

def movie(out):
 import sys
 sys.path.insert(0,os.environ.get('TOKYO_FFMPEG_DEPS','/pscratch/sd/d/dkololgi/graphweb_desi/outputs/tokyo_visuals_20261013_deps'))
 import imageio_ffmpeg
 plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
 style();d=np.load(out/'movie_galaxies.npz');xyz=d['xyz_Mpc'];p=d['P_WEB'].astype(float);cl=p.argmax(1);confidence=p.max(1)
 rgb=np.array([to_rgb(c) for c in CLASS_COLORS])[cl];rgba=np.column_stack([rgb,.2+.8*confidence]);center=(xyz.max(0)+xyz.min(0))/2;xyz=xyz-center
 f=plt.figure(figsize=(16,9),dpi=120);ax=f.add_axes([.02,.08,.75,.79],projection='3d',facecolor='black');ax.set_proj_type('ortho')
 spread=np.ptp(xyz,axis=0);ax.set_box_aspect(spread.copy());ax.set(xlim=(-spread[0]/2,spread[0]/2),ylim=(-spread[1]/2,spread[1]/2),zlim=(-spread[2]/2,spread[2]/2));ax.set_axis_off()
 points=ax.scatter(*xyz.T,s=2.2,c=rgba,depthshade=False,linewidths=0)
 f.text(.06,.92,'The inferred cosmic web in DESI',fontsize=28)
 f.text(.06,.872,'BGS galaxies · 0.20 ≤ z < 0.30 · observed redshift-space positions',fontsize=16,color='#C8CDD4')
 f.text(.785,.72,'Most probable\nenvironment',fontsize=20,linespacing=1.3)
 for j,(name,col) in enumerate(zip(CLASS_NAMES,CLASS_COLORS)):
  f.text(.795,.625-.07*j,'●',color=col,fontsize=25);f.text(.828,.625-.07*j,name,fontsize=19)
 f.text(.785,.29,'Brighter = more certain',fontsize=15)
 f.text(.785,.23,f'{len(xyz):,} displayed galaxies\nRandom galaxy subsample',fontsize=14,color='#C8CDD4',linespacing=1.5)
 f.text(.06,.055,'Provisional per-galaxy inference · camera rotation, not cosmic evolution',fontsize=15,color='#C8CDD4')
 f.text(.06,.025,'Frozen P12-A · real-DESI calibration unverified · distances in comoving Mpc',fontsize=13,color='#C8CDD4')
 # Box edges retain a scale reference without suggesting reconstructed field lines.
 lim=spread/2
 for axis in range(3):
  other=[j for j in range(3) if j!=axis]
  for a in [-1,1]:
   for b in [-1,1]:
    line=np.zeros((2,3));line[:,axis]=[-lim[axis],lim[axis]];line[:,other[0]]=a*lim[other[0]];line[:,other[1]]=b*lim[other[1]];ax.plot(*line.T,color='#708090',alpha=.3,lw=.5)
 f.text(.785,.14,f'Volume spans\n{spread[0]:.0f} × {spread[1]:.0f} × {spread[2]:.0f} Mpc',fontsize=14,color='#C8CDD4')
 frames=480;fps=24;writer=FFMpegWriter(fps=fps,codec='libx264',bitrate=7000,extra_args=['-pix_fmt','yuv420p','-threads','4','-movflags','+faststart'])
 with writer.saving(f,str(out/'DESI_environment_rotation.mp4'),dpi=120):
  for i in range(frames):
   t=i/fps;fade=np.clip((t-1)/2,0,1);c=rgba.copy();c[:,:3]=(1-fade)*.88+fade*rgb;points.set_facecolors(c)
   angle=-65+360*np.clip((t-3)/16,0,1);ax.view_init(elev=19,azim=angle)
   if i in [0,96,240,384]:f.savefig(out/f'movie_frame_{i:03d}.png',dpi=80)
   writer.grab_frame()
   if i%96==0:print('frame',i,'/',frames,flush=True)
 plt.close(f)
 previews=[Image.open(out/f'movie_frame_{i:03d}.png').convert('RGB') for i in [0,96,240,384]]
 contact=Image.new('RGB',(1280,720));
 for im,(x,y) in zip(previews,[(0,0),(640,0),(0,360),(640,360)]):contact.paste(im.resize((640,360)),(x,y))
 contact.save(out/'movie_contact_sheet.png')

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=Path,required=True);a.add_argument('--stage',choices=['prepare','figure','movie','render','all'],default='all');args=a.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 if args.stage in ['prepare','all']:prepare(args.out)
 if args.stage in ['figure','render','all']:figure(args.out)
 if args.stage in ['movie','render','all']:movie(args.out)
