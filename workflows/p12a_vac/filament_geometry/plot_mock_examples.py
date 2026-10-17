"""Matched observational-coordinate panels with explicit per-galaxy truth/predictions."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from prepare_parent import digest
from mock_truth import R
from compare_offslab import contrast2

COLORS=['#39A98F','#4089C0','#D65383','#C89417'];NAMES=['Void','Wall','Filament','Cluster'];MARKERS=['o','s','^','D']

def main():
 root=Path(sys.argv[1]);m=json.loads((root/'EXAMPLES.json').read_text());files=[]
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.12,'figure.facecolor':'white','savefig.facecolor':'white'})
 for row in m['examples']:
  key=row['key'];path=root/(key+'.npz');assert digest(path)==row['array_sha256'];d=np.load(path)
  xyz=d['delta_Mpc'];t=d['truth_class'];p=d['probability'];central=d['central'];known=np.isfinite(p).all(1);pred=np.full(len(p),-1);pred[known]=p[known].argmax(1)
  assert len(xyz)==row['display_rows'] and central.sum()==1
  slab=(abs(xyz[:,2])<10)&((xyz[:,:2]**2).sum(1)<=R*R)&~central
  np.testing.assert_allclose(contrast2(xyz[slab,:2]),row['projected_anisotropy'],atol=1e-10,rtol=0)
  np.testing.assert_array_equal(p[central][0],row['probability']);np.testing.assert_allclose(d['truth'][central][0],row['truth_eigenvalues'],rtol=0,atol=1e-6)
  fig=plt.figure(figsize=(14,10));gs=fig.add_gridspec(3,3,height_ratios=[1,1,.56],left=.07,right=.97,bottom=.105,top=.82,wspace=.22,hspace=.37)
  for view,(a,b,label) in enumerate([(0,1,'20 Mpc slab: projected view'),(0,2,'Full R sphere: orthogonal view')]):
   for col,title in enumerate(['Galaxy distribution','Native tidal truth at each galaxy','P12-A most probable class']):
    use=((abs(xyz[:,2])<10)&((xyz[:,:2]**2).sum(1)<=R*R)) if view==0 else (np.linalg.norm(xyz,axis=1)<=R)
    ax=fig.add_subplot(gs[view,col]);x=xyz[use,a];y=xyz[use,b]
    if col==0:ax.scatter(x,y,c='#767676',s=30,alpha=.75,edgecolors='white',linewidths=.3)
    else:
     values=t[use] if col==1 else pred[use]
     for j in range(4):
      k=values==j;alpha=np.ones(k.sum())*.9 if col==1 else .3+.7*p[use][k].max(1)
      rgba=np.tile(matplotlib.colors.to_rgba(COLORS[j]),(k.sum(),1));rgba[:,3]=alpha
      ax.scatter(x[k],y[k],c=rgba,s=38,marker=MARKERS[j],edgecolors='#333333',linewidths=.3)
     k=values<0;ax.scatter(x[k],y[k],c='#aaaaaa',s=28,marker='x',linewidths=.8)
    ax.scatter([0],[0],s=190,marker='*',facecolors='none',edgecolors='black',linewidths=1.2,zorder=10)
    ax.add_patch(Circle((0,0),R,fill=False,color='#555555',ls=':',lw=1))
    ax.set(xlim=(-1.1*R,1.1*R),ylim=(-1.1*R,1.1*R),aspect='equal',xlabel='Δx [Mpc]',ylabel=f'Δ{"y" if b==1 else "height"} [Mpc]')
    ax.set_xticks([-10,0,10]);ax.set_yticks([-10,0,10])
    if view==0:ax.set_title(title,fontsize=12,pad=9)
    if col==0:ax.text(.02,.98,label,transform=ax.transAxes,ha='left',va='top',fontsize=9)
    if view==1:
     ax.axhline(-10,color='#999999',ls='--',lw=.7);ax.axhline(10,color='#999999',ls='--',lw=.7)
  ax=fig.add_subplot(gs[2,0]);lam=np.array(row['truth_eigenvalues']);ax.bar([1,2,3],lam,color=['#777777' if v<=.2 else COLORS[int((lam>.2).sum())] for v in lam],width=.55)
  ax.axhline(.2,color='black',ls=':',label='Threshold 0.2');ax.set_xticks([1,2,3],['λ₁','λ₂','λ₃']);ax.set_title('Central galaxy: true eigenvalues',fontsize=11);ax.legend(fontsize=8)
  ax=fig.add_subplot(gs[2,1]);ax.bar(range(4),row['probability'],color=COLORS,width=.65);ax.set_xticks(range(4),NAMES,fontsize=9);ax.set(ylim=(0,1),title='Central galaxy: original probabilities')
  for j,v in enumerate(row['probability']):ax.text(j,v+.025,f'{v:.3f}',ha='center',fontsize=9)
  ax=fig.add_subplot(gs[2,2]);ax.axis('off');dominant=['Line','Sheet','Isotropic'][np.argmax(row['LPS_R'])]
  text=(f"Full 3D neighbourhood at R\nN = {row['neighbors_R']} (excluding centre)\nL / P / S = "+' / '.join(f'{v:.2f}' for v in row['LPS_R'])+f"\n{dominant}-dominated tracer covariance\n\nSame galaxies across columns in each row" )
  ax.text(0,.98,text,va='top',fontsize=10,linespacing=1.5)
  fig.suptitle(row['title'],fontsize=17,y=.98)
  fig.text(.5,.935,f"ph006 development mock · target {row['targetid']} · frozen R = 7 Mpc/h",ha='center',fontsize=11)
  handles=[Line2D([],[],marker=MARKERS[j],linestyle='',color=COLORS[j],label=NAMES[j],markersize=7) for j in range(4)]
  handles += [Line2D([],[],marker='x',linestyle='',color='#aaa',label='Unavailable',markersize=7),Line2D([],[],marker='*',linestyle='',markerfacecolor='none',color='black',label='Central galaxy',markersize=10)]
  fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.917),ncol=6,frameon=False,fontsize=10)
  fig.text(.5,.055,'Observed redshift-space coordinates. Top: 20 Mpc slab within projected radius R. Middle: full 3D sphere of radius R.\nDotted circles mark R; dashed lines mark slab boundaries. Prediction opacity increases with confidence; colours show argmax, not certainty.',ha='center',fontsize=9)
  fig.text(.5,.018,'Median parent ID within each predeclared category; illustrative cases, not prevalence estimates. Native host labels are not a displayed matter field.',ha='center',fontsize=9)
  for ext in ['png','pdf']:
   out=root/(key+'.'+ext);fig.savefig(out,dpi=170);files.append(out)
  plt.close(fig)
 validation={'pass':True,'examples':len(m['examples']),'same_rows_across_columns_per_view':True,'central_saved_probability_replay':True,'projected_anisotropy_replay':True,'neighbourhood_geometry_replay':'verified during extraction at R and 1.5R','plot_script_sha256':digest(__file__),'files':{p.name:digest(p) for p in files}}
 (root/'FIGURE_VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n');print(json.dumps(validation),flush=True)
if __name__=='__main__':main()
