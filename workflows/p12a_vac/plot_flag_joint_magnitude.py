"""Shared-scale population diagrams and explicitly conditional retention plots."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,TwoSlopeNorm

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root
 d=np.load(r/'HISTOGRAMS.npz');z=d['z_edges'];m=d['r_edges'];area=d['area_deg2'];real=d['loa'];mock=d['mock_observed'];mid=(z[1:]+z[:-1])/2
 plt.rcParams.update({'font.size':10,'axes.titlesize':11,'figure.titlesize':14})
 fig,ax=plt.subplots(2,3,figsize=(14,8),sharex=True,sharey=True,layout='constrained')
 dens=np.stack([real/area[:,None,None],mock/area[:,None,None]])
 norm=LogNorm(vmin=.002,vmax=float(dens.max()))
 for c,name in enumerate(['SGC','NGC']):
  for j,label in enumerate(['Loa selected galaxies','Loa-processed ph006 mock']):
   h=dens[j,c];im=ax[c,j].pcolormesh(z,m,np.ma.masked_where(h.T<=0,h.T),norm=norm,cmap='viridis',rasterized=True)
   ax[c,j].set_title(name+' — '+label)
  ratio=np.divide(real[c],mock[c],out=np.full_like(real[c],np.nan,dtype=float),where=mock[c]>0)
  ratio[(real[c]<20)|(mock[c]<20)]=np.nan
  ri=ax[c,2].pcolormesh(z,m,ratio.T,cmap='coolwarm',norm=TwoSlopeNorm(vmin=0,vcenter=1,vmax=3),rasterized=True)
  ax[c,2].set_title(name+' — Loa / mock counts')
  for j in range(3):
   ax[c,j].axhline(19.5,color='black' if j==2 else 'white',ls='--',lw=.9)
   if c==1 and j==0:ax[c,j].axhline(19.54,color='black' if j==2 else 'white',ls=':',lw=.8)
   ax[c,j].set_ylim(20.2,12);ax[c,j].set_xlim(.15,.55);ax[c,j].set_xlabel('Redshift')
  ax[c,0].set_ylabel('Apparent r magnitude')
 fig.colorbar(im,ax=ax[:,:2],label='Galaxies / deg² per Δz=0.01, Δr=0.1 bin',shrink=.85)
 fig.colorbar(ri,ax=ax[:,2],label='Count ratio (cells with ≥20 in each sample)',extend='max',shrink=.85)
 fig.suptitle('Magnitude–redshift populations on identical angular support\nLoa: extinction-corrected Legacy r; mock: stored R_MAG_APP (passband equivalence unresolved)')
 fig.savefig(r/'magnitude_redshift.png',dpi=160);fig.savefig(r/'magnitude_redshift.pdf');plt.close(fig)
 fig,ax=plt.subplots(2,3,figsize=(14,8),layout='constrained')
 allrat=d['combinations']/np.maximum(mock.sum(-1)[...,None],1)
 for c,name in enumerate(['SGC','NGC']):
  for j,(num,den,title) in enumerate([(mock[c],d['mock_parent'][c],'Mock successful / forFA parent'),(d['loa_strict'][c],real[c],'Loa all eight cuts / baseline')]):
   frac=np.divide(num,den,out=np.full_like(num,np.nan,dtype=float),where=den>=20)
   im=ax[c,j].pcolormesh(z,m,frac.T,vmin=0,vmax=1,cmap='cividis',rasterized=True);ax[c,j].set_ylim(20.2,12);ax[c,j].set_xlabel('Redshift');ax[c,j].set_title(name+' — '+title)
  ax[c,0].set_ylabel('Apparent r magnitude')
  ax[c,2].fill_between(mid,allrat[c].min(-1),allrat[c].max(-1),alpha=.25,color='tab:blue',label='Envelope of 256 combinations')
  ax[c,2].plot(mid,allrat[c,:,0],label='Baseline',color='tab:blue')
  ax[c,2].plot(mid,allrat[c,:,255],label='All eight cuts',color='tab:orange',ls='--')
  ax[c,2].axhline(1,color='black',lw=.8);ax[c,2].set_xlabel('Redshift');ax[c,2].set_ylabel('Loa / observed mock counts');ax[c,2].set_title(name+' — joint-cut effect on N(z)');ax[c,2].set_ylim(0,allrat.max()*1.05);ax[c,2].grid(alpha=.2)
 ax[0,2].legend(fontsize=8)
 fig.colorbar(im,ax=ax[:,:2],label='Conditional retention fraction; denominator ≥20',shrink=.85)
 fig.suptitle('Selection retention and all joint cut combinations\nMock retention includes processing/assignment losses; Loa retention is not absolute survey completeness')
 fig.savefig(r/'retention_and_joint_cuts.png',dpi=160);fig.savefig(r/'retention_and_joint_cuts.pdf');plt.close(fig)
if __name__=='__main__':main()
