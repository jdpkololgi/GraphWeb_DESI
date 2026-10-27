#!/usr/bin/env python3
"""Bounded original-slab row audit and honest filament presentation. No inference."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,hashlib,shutil
import numpy as np, fitsio
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
ROOT=Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010')
OUT=ROOT/'filament_and_local_geometry_20261017'
WHITE,GREY='#F5F5F5','#B9C2D0'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def export(f,n):
 for ext in ('png','pdf','svg'):f.savefig(OUT/(n+'.'+ext),dpi=240 if ext=='png' else 120,facecolor='black')
 plt.close(f)
 assert b'/FontFile2' in (OUT/(n+'.pdf')).read_bytes()
def main():
 for d in ('data','scripts'):(OUT/d).mkdir(parents=True,exist_ok=True)
 source=ROOT/'a06_desi';m=json.loads((source/'manifest.json').read_text());path=source/'plotted_galaxies.npz'
 assert sha(path)==m['output_hashes']['plotted_galaxies.npz']
 with np.load(path) as z:d={k:z[k] for k in z.files}
 with fitsio.FITS(m['catalogue_path']) as f:a=f[1].read(rows=d['catalogue_row'],columns=['TARGETID','RA','DEC','Z','P_WEB','SUPPORTED','QUALITY','BOUNDARY_MPC','EIGENVALUE_MEAN'])
 for k in ['TARGETID','RA','DEC','Z','SUPPORTED','QUALITY']:np.testing.assert_array_equal(a[k],d[k])
 p=a['P_WEB'].astype(float);np.testing.assert_array_equal(p[:,2],d['P_filament']);np.testing.assert_allclose(p.sum(1),1,atol=1e-6)
 x,y=d['x_Mpc'],d['y_Mpc'];xy=np.column_stack((x,y));r=7/.6766
 # Diagnostic on displayed projected points only, not a 3D filament truth label.
 tree=cKDTree(xy);anis=np.full(len(x),np.nan);count=np.zeros(len(x),int)
 for i,near in enumerate(tree.query_ball_point(xy,r)):
  near=[j for j in near if j!=i];count[i]=len(near)
  if len(near)>=8:
   eig=np.linalg.eigvalsh(np.cov(xy[near].T));anis[i]=(eig[1]-eig[0])/eig.sum()
 interior=(x>-180+r)&(x<180-r)&(y>800+r)&(y<1250-r)
 valid=interior&np.isfinite(anis);threshold=float(np.quantile(anis[valid],.8));linear=valid&(anis>=threshold)
 low=p[:,2]<.2
 def group(mask):return dict(n=int(mask.sum()),mean_probabilities=p[mask].mean(0).tolist(),argmax_counts=np.bincount(p[mask].argmax(1),minlength=4).tolist(),low_filament_count=int((mask&low).sum()))
 metrics=dict(all=group(np.ones(len(x),bool)),low_filament=group(low),projected_line_like=group(linear),projected_line_like_low_filament=group(linear&low),
   local_shape=dict(radius_mpc=r,min_neighbors=8,valid_rows=int(valid.sum()),top20pct_cut=threshold,definition='Covariance eigenvalue contrast of neighbouring displayed projected points; excludes self and plot-edge targets; not independent truth or full 3D environment'),
   quality_values=np.unique(a['QUALITY']).tolist(),boundary_min_mpc=float(a['BOUNDARY_MPC'].min()),
   integrity='All original rows, coordinates, IDs, support/quality and P(filament) exactly match frozen VAC',
   source=dict(path=m['catalogue_path'],sha256_inherited=m['catalogue_sha256'],plotted_rows_sha256=sha(path)),
   smoothing_radius_mpc=r,slice_depth_mpc=20,interpretation='T-web local eigenvalue classes, not a geometric filament skeleton. Shape diagnostic is descriptive and dependent; cannot identify the cause of disagreement or validate real-data calibration.')
 write(OUT/'data/filament_audit.json',metrics)
 np.savez_compressed(OUT/'data/filament_audit_rows.npz',**d,P_WEB=p,eigenvalue_mean=a['EIGENVALUE_MEAN'],projected_anisotropy=anis,projected_neighbor_count=count,projected_shape_valid=valid,projected_line_like=linear)
 plt.rcParams.update({'font.family':'DejaVu Sans','text.color':WHITE,'axes.labelcolor':GREY,'xtick.color':GREY,'ytick.color':GREY,'axes.edgecolor':'#626C7A','pdf.fonttype':42,'svg.fonttype':'none'})
 f=plt.figure(figsize=(16,9),dpi=120,facecolor='black');names=['Void','Wall','Filament','Cluster']
 for k in range(4):
  ax=f.add_axes([.055+.235*k,.22,.205,.62],facecolor='black')
  im=ax.scatter(x,y,c=p[:,k],s=3,cmap='viridis',vmin=0,vmax=1,lw=0)
  ax.set(xlim=(-180,180),ylim=(800,1250),aspect='equal',xticks=[-150,0,150],yticks=[850,1000,1150])
  if k:ax.tick_params(labelleft=False)
  else:ax.set_ylabel('Radial-plane coordinate [Mpc]',fontsize=14)
  ax.set_title('P('+names[k].lower()+')',fontsize=20,pad=12)
  ax.set_xlabel('Transverse [Mpc]',fontsize=13);ax.tick_params(labelsize=12)
  for sp in ('top','right'):ax.spines[sp].set_visible(False)
 cax=f.add_axes([.30,.12,.4,.018]);cb=f.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,.25,.5,.75,1]);cb.set_label('Saved class probability',fontsize=14)
 f.text(.5,.94,'Where does the probability go?',ha='center',fontsize=27)
 f.text(.5,.045,'Identical 4,459 galaxies in every panel · 20 Mpc slice · fixed 0–1 scales · provisional P12-A',ha='center',fontsize=14,color=GREY)
 export(f,'four_environment_probabilities')
 # Geometric centre crop, no optimization for P(filament) or attractive morphology.
 crop=(abs(x)<100)&(y>=925)&(y<1125);high=crop&(p[:,2]>=.5)
 f=plt.figure(figsize=(16,9),dpi=120,facecolor='black');ax=f.add_axes([.07,.14,.56,.76],facecolor='black')
 ax.scatter(x[crop],y[crop],s=8,c='#858D99',alpha=.65,lw=0)
 cmap=LinearSegmentedColormap.from_list('filament',['#E85D9E','#FFF0A5'])
 im=ax.scatter(x[high],y[high],s=13,c=p[high,2],cmap=cmap,vmin=.5,vmax=1,lw=0)
 ax.set(xlim=(-100,100),ylim=(925,1125),aspect='equal',xticks=[-100,-50,0,50,100],yticks=[950,1000,1050,1100]);ax.tick_params(labelsize=13)
 ax.set_xlabel('Transverse coordinate [Mpc]',fontsize=16);ax.set_ylabel('Radial-plane coordinate [Mpc]',fontsize=16)
 ax.add_patch(Circle((-80,1105),r,fill=False,edgecolor=WHITE,lw=1.2));ax.text(-63,1105,'R = 7 h⁻¹ Mpc',fontsize=12,color=WHITE,va='center')
 f.text(.655,.83,'Filament-probable galaxies',fontsize=23)
 f.text(.655,.72,f'{int(high.sum()):,} with P(filament) ≥ 0.5',fontsize=18,color='#F1A8CC')
 f.text(.655,.64,f'All {int(crop.sum()):,} galaxies retained as context',fontsize=16,color=GREY)
 cax=f.add_axes([.66,.55,.285,.022]);cb=f.colorbar(im,cax=cax,orientation='horizontal',ticks=[.5,.75,1]);cb.set_label('P(filament)',fontsize=15)
 f.text(.655,.405,'Local tidal classification',fontsize=19)
 f.text(.655,.345,'R = 7 h⁻¹ Mpc · λₜₕ = 0.2',fontsize=16,color=GREY)
 f.text(.655,.255,'Not an extracted filament network',fontsize=16,color=GREY)
 f.text(.655,.195,'Central crop of the original slab',fontsize=15,color=GREY)
 f.text(.5,.045,'Observed redshift-space positions · 20 Mpc depth · frozen provisional catalogue · no field reconstruction',ha='center',fontsize=13,color=GREY)
 export(f,'filament_probability_highlight')
 metrics['new_plot']=dict(crop_x=[-100,100],crop_y=[925,1125],crop_rule='Geometric central 200x200 Mpc of original bounds, no optimized region search',galaxies=int(crop.sum()),highlighted=int(high.sum()),threshold=.5,background_retained=True)
 write(OUT/'data/filament_audit.json',metrics);shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
 print(json.dumps(metrics))
if __name__=='__main__':main()
