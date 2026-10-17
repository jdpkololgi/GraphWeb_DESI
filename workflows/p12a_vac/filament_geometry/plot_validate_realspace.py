"""Independent numeric checks and paired matter/galaxy/RSD visual comparisons."""
import json,sys
from pathlib import Path
import numpy as np
from compare_offslab import shape
from prepare_parent import digest
from mock_truth import R


def svd_shape(x,w=None):
 if w is None:w=np.ones(len(x))
 w=w/w.sum();center=np.sum(x*w[:,None],axis=0);s=np.linalg.svd((x-center)*np.sqrt(w[:,None]),compute_uv=False)**2
 return np.array([(s[0]-s[1])/s[0],(s[1]-s[2])/s[0],s[2]/s[0]])


def main():
 root=Path(sys.argv[1]);report=json.loads((root/'COMPARISON.json').read_text());checks=0
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import Normalize
 from matplotlib.patches import Circle
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':False})
 colors=['#bc5b22','#367fb1','#888888'];files=[]
 for row in report['regions']:
  key=row['key'];path=root/(key+'_comparison.npz');assert digest(path)==row['arrays_sha256'];d=np.load(path);real=d['real_delta_Mpc'];rs=d['rs_delta_Mpc'];exact=d['exact_observed_rs_delta_Mpc'];obs=d['observed'];inreal=d['inreal'];inrs=d['inrs'];central=d['central'];ids=d['targetid']
  assert central.sum()==1 and inrs.sum()==row['observed_count_target'];checks+=2
  row['raw_vs_observed_rs_max_difference_Mpc']=float(np.max(np.linalg.norm(rs[obs]-exact[obs],axis=1)))
  masks={'all_available_galaxies_real':(real,inreal),'bright_parent_real':(real,inreal&(ids>0)),'observed_fixed_real_ball':(real,inreal&obs),'observed_fixed_real_ids_in_rs':(exact,inreal&obs),'observed_fixed_rs_ball':(exact,inrs),'observed_fixed_rs_ids_in_real':(real,inrs)}
  for name,(x,k) in masks.items():
   assert k.sum()==row['galaxies'][name]['n'];np.testing.assert_allclose(svd_shape(x[k]),row['galaxies'][name]['LPS'],atol=1e-12);checks+=2
  for name,k in [('all_available_galaxies',inreal),('bright_parent',inreal&(ids>0))]:
   if row['thinning'][name]['available']:
    ix=d[name+'_draw_indices'];assert ix.shape==(512,int(inrs.sum()));assert np.all(np.sort(ix,axis=1)[:,1:]!=np.sort(ix,axis=1)[:,:-1]);checks+=2
    for j in range(512):
     np.testing.assert_allclose(svd_shape(real[k][ix[j]]),d[name+'_thin_real_LPS'][j],atol=1e-12)
     np.testing.assert_allclose(svd_shape(rs[k][ix[j]]),d[name+'_thin_rs_LPS'][j],atol=1e-12);checks+=2
    assert sum(row['thinning'][name]['real']['dominant_counts'])==512;checks+=1
  axis=d['grid_axis_Mpc'];xx=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1);inside=np.linalg.norm(xx,axis=3)<=R;mass=d['matter_counts'];assert np.isfinite(mass).all() and (mass>=0).all();checks+=1
  np.testing.assert_allclose(svd_shape(xx[inside],mass[inside]),row['matter']['total_mass_LPS'],atol=2e-7)
  np.testing.assert_allclose(svd_shape(xx[inside],np.maximum(mass[inside]-d['cosmic_mean'],0)),row['matter']['positive_excess_LPS'],atol=2e-7);checks+=2
  for keya in ['matter_thin_LPS','all_available_galaxies_thin_real_LPS','all_available_galaxies_thin_rs_LPS']:
   a=d[keya];assert a.shape==(512,3);np.testing.assert_allclose(a.sum(1),1,atol=1e-12);checks+=2
  np.testing.assert_allclose(svd_shape(xx[inside]),[0,0,1],atol=1e-8);checks+=1
  thin_index=d['all_available_galaxies_draw_indices'][0];example=real[inreal][thin_index]
  samples=[None,real[inreal],example,real[inrs],exact[inrs]]
  specs=[row['matter']['total_mass_LPS'],row['galaxies']['all_available_galaxies_real']['LPS'],row['thinning']['all_available_galaxies']['real']['median_LPS'],row['galaxies']['observed_fixed_rs_ids_in_real']['LPS'],row['galaxies']['observed_fixed_rs_ball']['LPS']]
  titles=['Dense matter\nreal space',f"All available galaxies\nreal space · N={inreal.sum()}",f"Random thinning\nreal space · N={inrs.sum()}",f"Observed sample\nreal space · N={inrs.sum()}",f"Same observed IDs\nredshift space · N={inrs.sum()}"]
  fig=plt.figure(figsize=(19,11));gs=fig.add_gridspec(3,5,height_ratios=[1,1,.55],left=.055,right=.98,bottom=.18,top=.86,wspace=.22,hspace=.40)
  for view,(a,b,integrate) in enumerate([(0,1,2),(0,2,1)]):
   # Column density through the same real-space R sphere, with fixed cosmic
   # mean normalization; values below display floor remain in quantitative data.
   density_projection=np.sum(np.where(inside,mass,0),axis=integrate)*(axis[1]-axis[0])/(float(d['cosmic_mean'])*2*R)
   z=np.log10(np.maximum(density_projection,1e-3))
   for col in range(5):
    ax=fig.add_subplot(gs[view,col])
    if col==0:im=ax.imshow(z.T,origin='lower',extent=[axis[0],axis[-1],axis[0],axis[-1]],cmap='magma',norm=Normalize(-.5,1.5),interpolation='nearest')
    else:
     x=samples[col];ax.scatter(x[:,a],x[:,b],s=23,color='#285c81' if col<4 else '#9c4b3f',alpha=.8,edgecolor='white',linewidth=.25)
    ax.add_patch(Circle((0,0),R,fill=False,color='white' if col==0 else '#888',ls=':',lw=1))
    ax.scatter([0],[0],s=95,marker='*',facecolors='none',edgecolors='white' if col==0 else 'black',linewidths=1,zorder=10)
    limit=1.5*R
    if col>0:assert np.all(abs(samples[col][:,[a,b]])<=limit)
    ax.set(xlim=(-limit,limit),ylim=(-limit,limit),aspect='equal',xlabel='Δx [Mpc]');ax.set_xticks([-10,0,10]);ax.set_yticks([-10,0,10])
    if col==0:ax.set_ylabel('Δy [Mpc]' if view==0 else 'Δheight [Mpc]')
    if view==0:ax.set_title(titles[col],fontsize=11,pad=8)
  for col in range(5):
   ax=fig.add_subplot(gs[2,col]);v=np.asarray(specs[col]);ax.bar(range(3),v,color=colors,width=.6)
   if col==2:
    q=np.asarray(row['thinning']['all_available_galaxies']['real']['q025_q975_LPS']);ax.errorbar(range(3),v,yerr=[v-q[0],q[1]-v],fmt='none',color='black',capsize=3,lw=1)
   for j,value in enumerate(v):ax.text(j,min(value+.025,.96),f'{value:.2f}',ha='center',fontsize=9)
   ax.set_xticks(range(3),['Linearity','Planarity','Sphericity'],rotation=25,fontsize=8);ax.set_ylim(0,1)
   ax.set_title('512-draw median + 95% range' if col==2 else 'Full 3D covariance',fontsize=9)
  fig.suptitle(f"Same mock region · TARGETID {row['targetid']} · R = 7/0.6766 Mpc",fontsize=16,y=.97)
  fig.text(.5,.925,'Matter → galaxies → random thinning → actual observation selection → redshift distortions',ha='center',fontsize=12)
  cax=fig.add_axes([.065,.065,.16,.012]);fig.colorbar(im,cax=cax,orientation='horizontal',label='log₁₀[Σ / (2R × cosmic mean density)]');cax.tick_params(labelsize=8);cax.xaxis.label.set_size(8)
  fig.text(.62,.065,'Fixed first thinning draw shown; bars summarize 512 draws. Thin and observed samples need not contain the same IDs.\nLast two columns use identical IDs selected in the observed redshift-space R sphere; real positions may lie outside R.\n10% particle TSC matter grid, 0.977 Mpc/h resolution. These are tracer geometry statistics, not tidal classes.',ha='center',fontsize=9)
  for ext in ['png','pdf']:
   f=root/(row['key']+'_realspace.'+ext);fig.savefig(f,dpi=170);files.append(f)
  plt.close(fig)
 validation={'pass':True,'checks':checks,'script_sha256':digest(__file__),'comparison_sha256':digest(root/'COMPARISON.json'),'figure_sha256':{p.name:digest(p) for p in files},'raw_vs_observed_rs_max_difference_Mpc':{a['key']:a['raw_vs_observed_rs_max_difference_Mpc'] for a in report['regions']}}
 (root/'VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n');print(json.dumps(validation),flush=True)
if __name__=='__main__':main()
