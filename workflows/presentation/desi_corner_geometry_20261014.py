#!/usr/bin/env python3
"""Explain ordered-eigenvalue corner geometry from existing full-range histograms."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
BASE=Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010')
SRC=BASE/'desi_visuals_20261013/corner_sweep_arrays.npz'
OUT=BASE/'density_checks_20261014'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
a=np.load(SRC);e=a['full_edges'];c=(e[:-1]+e[1:])/2;x,y=np.meshgrid(c,c,indexing='ij');center=(x+y)/2;gap=y-x
h=a['full_pair_histograms'];names=['Abacus truth','Abacus posterior','DESI posterior'];col=['#CCD2DA','#5B8FF9','#F04D9A'];stats={}
plt.rcParams.update({'font.size':12,'figure.facecolor':'#080C13','axes.facecolor':'#080C13','text.color':'#EDF1F5','axes.labelcolor':'#EDF1F5','xtick.color':'#BBC5D2','ytick.color':'#BBC5D2','axes.edgecolor':'#8995A5','pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'#080C13'})
f,ax=plt.subplots(2,3,figsize=(15,10));f.subplots_adjust(left=.07,right=.97,bottom=.16,top=.84,wspace=.3,hspace=.40)
def contour(ax,z,xe,ye,color):
 z=gaussian_filter(z,1);s=np.sort(z.ravel())[::-1];q=np.cumsum(s);lev=np.unique([s[np.searchsorted(q,t*q[-1])] for t in [.9,.5]])
 ax.contour(xe,ye,z.T,levels=lev,colors=color,linewidths=1.8)
for m,n in enumerate(names):
 stats[n]={}
 for p,pair in enumerate(['lambda1_lambda2','lambda1_lambda3','lambda2_lambda3']):
  w=h[m,p]/h[m,p].sum()
  def var(v):mu=(w*v).sum();return (w*(v-mu)**2).sum()
  vx=var(x);vy=var(y);cov=(w*(x-(w*x).sum())*(y-(w*y).sum())).sum()
  stats[n][pair]={'binned_pearson':float(cov/np.sqrt(vx*vy)),'parallel_to_perpendicular_variance_ratio':float(var((x+y)/np.sqrt(2))/var((y-x)/np.sqrt(2))),'mass_below_ordering_boundary':float(w[y<x].sum())}
 z=h[m,1];contour(ax[0,m],z,c,c,col[m]);ax[0,m].plot([-.6,1.8],[-.6,1.8],':',color='#9AA4B4');ax[0,m].set(xlim=(-.6,1.8),ylim=(-.6,1.8),xlabel='λ₁',ylabel='λ₃',title=n);ax[0,m].set_aspect('equal')
 ax[0,m].text(.05,.93,'λ₃ ≥ λ₁: ordering boundary',transform=ax[0,m].transAxes,fontsize=10)
 ce=np.linspace(-.6,1.8,100);ge=np.linspace(0,2.4,100);hh=np.histogram2d(center.ravel(),gap.ravel(),bins=[ce,ge],weights=z.ravel())[0]
 contour(ax[1,m],hh,(ce[:-1]+ce[1:])/2,(ge[:-1]+ge[1:])/2,col[m]);ax[1,m].set(xlim=(-.6,1.8),ylim=(0,2.4),xlabel='Pair midpoint (λ₁ + λ₃) / 2',ylabel='Eigenvalue gap λ₃ − λ₁')
 for aa in ax[:,m]:aa.grid(alpha=.12);aa.spines[['top','right']].set_visible(False)
f.suptitle('Why the eigenvalue corner contours follow the diagonal',fontsize=24,y=.95)
f.text(.5,.90,'Ordering restricts the allowed half-plane; a common eigenvalue shift runs parallel to λ₃ = λ₁.',ha='center',fontsize=13)
f.text(.5,.085,'Bottom: the same binned distribution in midpoint–gap coordinates · gap measures anisotropy, midpoint includes trace and shape',ha='center',fontsize=12)
f.text(.5,.043,'Existing saved populations; no new inference · 50/90% display contours · binned transformation · native selections differ',ha='center',fontsize=11,color='#BBC5D2')
for ext in ['png','pdf','svg']:f.savefig(OUT/f'eigenvalue_geometry.{ext}',dpi=170,bbox_inches='tight')
plt.close(f);im=Image.open(OUT/'eigenvalue_geometry.png');im.thumbnail((1400,1000));im.save(OUT/'eigenvalue_geometry_preview.png')
(OUT/'eigenvalue_geometry_statistics.json').write_text(json.dumps({'source':str(SRC),'source_sha256':sha(SRC),'script_sha256':sha(__file__),'method':'Moments approximate from full-range 200x200 bin centres; no new posterior draws. Ratio is variance along (x+y)/sqrt2 divided by variance along (y-x)/sqrt2. Display smoothing may spill over ordering boundary. Midpoint-gap contours rebin existing histograms.','results':stats},indent=2)+'\n')
print(json.dumps(stats,indent=2))
