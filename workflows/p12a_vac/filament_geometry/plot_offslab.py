"""Plot registered geometry diagnostics; no model inference or catalogue edits."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from prepare_parent import digest
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from shared.plot_style import apply_style

def main():
 p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--parent',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 r=json.loads((a.results/'RESULTS.json').read_text());assert digest(a.results/'geometry.npz')==r['provenance']['geometry_sha256'];g=dict(np.load(a.results/'geometry.npz'));apply_style()
 cohort=g['cohort_low_fil_wall'];mask=cohort&g['eligible_1']&(g['nfull_1']>=8)&(g['nslab_1']>=8);f,axs=plt.subplots(2,2,figsize=(13,10),layout='constrained')
 for ax,key,title in zip(axs[0],['slab_1','full_1'],['Neighbours inside original slab','Neighbours in full 3D sphere']):
  v=g[key][mask];ax.scatter(v[:,0],v[:,1],s=13,alpha=.7,c='#4E84F7');ax.plot([0,1],[0,1],color='white',lw=.8,ls='--');ax.set(xlim=(0,1),ylim=(0,1),xlabel='Linearity L',ylabel='Planarity P',title=title);ax.set_aspect('equal')
 colors=['#F5C144','#4E84F7','#A1FCDD']
 for key,label,col in zip(['slab_1','full_1','matched_1'],['Slab','Full sphere','Full, matched neighbour count'],colors):
  v=g[key][mask];axs[1,0].hist(v[:,1]-v[:,0],bins=np.linspace(-1,1,26),histtype='step',density=True,label=label,color=col,lw=1.5)
 axs[1,0].set(xlabel='Planarity − linearity',ylabel='Density');axs[1,0].legend(fontsize=9)
 for name,col in [('low_fil_wall','#4E84F7'),('low_fil_void','#A1FCDD'),('line_fil_argmax','#EB336F')]:
  rows=[v for v in r['rows'] if v['cohort']==name and v['min_neighbors']==8];x=[v['radius_Mpc'] for v in rows];y=[v['full']['sheet_dominant']['mean'] for v in rows];ci=np.array([v['full']['sheet_dominant']['ci95'] for v in rows]);axs[1,1].errorbar(x,y,yerr=np.array([np.array(y)-ci[:,0],ci[:,1]-np.array(y)]),label=name,fmt='o-',color=col)
 axs[1,1].set(xlabel='Radius [Mpc]',ylabel='Fraction with P > L and P > S',ylim=(0,1));axs[1,1].legend(fontsize=9)
 f.suptitle(f'Off-slab tracer geometry: {mask.sum()} paired low-P(filament), wall-argmax targets\nCovariance shapes are not tidal truth; bars: conditional 95% spatial-block intervals',fontsize=13)
 for ext in ['png','pdf']:f.savefig(a.out/f'01_geometry_comparison.{ext}',dpi=170,bbox_inches='tight')
 plt.close(f)
 # Original orientation, offset0. Each series uses same all-depth eligible targets.
 f,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
 for name,col in [('low_fil_wall','#4E84F7'),('low_fil_void','#A1FCDD'),('line_fil_argmax','#EB336F')]:
  z=[v for v in r['views'] if v['cohort']==name and v['angle']==0 and v['offset']==0];axs[0].plot([v['depth'] for v in z],[v['line_like']/v['valid'] if v['valid'] else np.nan for v in z],'o-',label=name,color=col)
  common=g['cohort_'+name]&np.logical_and.reduce([np.isfinite(g[f'shape_a0_o0_d{dep}']) for dep in [20,40,80]]);axs[1].plot([20,40,80],[np.mean(g[f'shape_a0_o0_d{dep}'][common]>=.5847730874520835) for dep in [20,40,80]],'o-',label=f'{name} (n={common.sum()})',color=col)
 for ax,title in zip(axs,['Changing eligible targets','Same targets at 20 / 40 / 80 Mpc']):ax.set(xlabel='Full slab depth [Mpc]',ylabel='Fraction above fixed projected-line threshold',ylim=(0,1),title=title);ax.legend(fontsize=8)
 f.suptitle('Thickness changes both tracer eligibility and projected morphology',fontsize=14)
 for ext in ['png','pdf']:f.savefig(a.out/f'02_thickness.{ext}',dpi=170,bbox_inches='tight')
 plt.close(f)
 # Illustrations selected by TARGETID rank, without optimizing geometry or crops.
 d=dict(np.load(a.parent/'parent.npz'));assert digest(a.parent/'parent.npz')==r['provenance']['parent_sha256'];xyz=d['xyz_Mpc'];tree=cKDTree(xyz);rad=7/.6766
 pool=np.flatnonzero(cohort&g['eligible_1']&(g['nfull_1']>=8));pool=pool[np.argsort(g['TARGETID'][pool])];pick=pool[np.floor((len(pool)-1)*np.array([.25,.5,.75])).astype(int)];f,axs=plt.subplots(3,3,figsize=(13,12),layout='constrained');examples=[]
 for row,i in enumerate(pick):
  j=d['original_target_parent_index'][i];q=np.array(tree.query_ball_point(xyz[j],rad));q=q[q!=j];delta=xyz[q]-xyz[j];ins=abs(xyz[q,2])<10
  for ax,(u,v),title in zip(axs[row],[(0,1),(0,2),(1,2)],['x–y (original viewing plane)','x–height','y–height']):
   ax.scatter(delta[ins,u],delta[ins,v],s=16,color='#F5C144',label='Inside original slab');ax.scatter(delta[~ins,u],delta[~ins,v],s=16,color='#A1FCDD',marker='x',label='Off slab');ax.scatter([0],[0],marker='*',color='white',s=100)
   ax.set(xlim=(-rad,rad),ylim=(-rad,rad),xlabel=['Δx','Δy','Δheight'][u]+' [Mpc]',ylabel=['Δx','Δy','Δheight'][v]+' [Mpc]');ax.set_aspect('equal');ax.set_title(title,fontsize=11)
  prob=g['P_WEB'][i];axs[row,0].text(.01,1.15,f'TARGETID {g["TARGETID"][i]} · P(wall)={prob[1]:.2f}, P(fil)={prob[2]:.2f}',transform=axs[row,0].transAxes,fontsize=10)
  examples.append({'TARGETID':int(g['TARGETID'][i]),'target_rank_quantile':[.25,.5,.75][row],'neighbors':len(q),'offslab':int((~ins).sum()),'LPS':g['full_1'][i].tolist()})
 axs[0,2].legend(fontsize=8,loc='upper right');f.suptitle('Same neighbours viewed in three planes\nExamples chosen by TARGETID rank; no favourable morphology selection',fontsize=14)
 for ext in ['png','pdf']:f.savefig(a.out/f'03_neighbour_views.{ext}',dpi=170,bbox_inches='tight')
 plt.close(f)
 (a.out/'EXAMPLES.json').write_text(json.dumps(examples,indent=2)+'\n')
 (a.out/'ARTIFACT_SHA256.json').write_text(json.dumps({p.name:digest(p) for p in a.out.iterdir() if p.suffix in ['.png','.pdf','.json'] and p.name!='ARTIFACT_SHA256.json'},indent=2)+'\n')
if __name__=='__main__':main()
