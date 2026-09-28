"""Plot preserved full-catalogue counts; single-phase diagnostic, no error bands."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from shared.plot_style import apply_style,finalize_axes
apply_style()
x=json.load(open(R/'docs/evidence/p12a_loa_completeness_20260928/RESULTS.json'))
p=R/'docs/evidence/p12a_alignment_execution_20260926/parent'
z=np.linspace(.155,.545,40);v1=np.load(p/'v1.npz')['z'];old=np.load(p/'v0.1.npz')['z']
obs=np.array(x['fine_counts']['uniform_r19p5_observed']);correct=np.array(x['fine_counts']['uniform_r19p5_assignment_zfail'])
f,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey='row',layout='constrained')
for i,cap in enumerate(['SGC','NGC']):
 for y,label,ls in [(old,'Old raw parent','--'),(v1,'v1 raw parent','-'),(obs,'Loa observed',':'),(correct,'Loa assignment + zfail corrected','-.')]:axes[0,i].plot(z,y[i],label=label,ls=ls)
 axes[0,i].set(title=cap,ylabel='Count per dz = 0.01',yscale='log')
 for y,label,ls in [(old,'Old / corrected Loa','--'),(v1,'v1 / corrected Loa','-')]:axes[1,i].plot(z,y[i]/correct[i],label=label,ls=ls)
 axes[1,i].axhline(1,color='white',ls=':',lw=1);axes[1,i].set(xlabel='Redshift',ylabel='Raw parent / corrected Loa',ylim=(0,2.5))
 for a in axes[:,i]:finalize_axes(a,a.get_title(),a.get_xlabel(),a.get_ylabel())
axes[0,0].legend(fontsize=9);axes[1,0].legend(fontsize=9)
f.suptitle('Common sky, 12 ≤ r < 19.5; one ph000 parent\nPixel support only: exact imaging/veto parity remains under test',fontsize=13)
out=R/'docs/figures/p12a_loa_completeness_20260928';out.mkdir(exist_ok=True)
for ext in ['png','pdf']:f.savefig(out/('counts.'+ext),dpi=160)
