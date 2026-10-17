"""Independent saved-array reconciliation and compact evidence figure."""
import json,sys
from pathlib import Path
import numpy as np
from prepare_parent import digest

def main():
 root=Path(sys.argv[1]); r=json.loads((root/'RESULTS.json').read_text());d=np.load(root/'geometry.npz')
 assert digest(root/'geometry.npz')==r['geometry_sha256']
 p=d['probability'];t=(d['truth_eigenvalues']>.2).sum(1);w=d['natural_weight'];n=len(t)
 assert len(np.unique(d['parent_node_id']))==n
 assert np.isfinite(p).all() and (p>=0).all() and (p<=1).all()
 assert np.isfinite(w).all() and (w>0).all()
 assert (np.diff(d['truth_eigenvalues'],axis=1)>=0).all()
 np.testing.assert_allclose(p.sum(1),1)
 checks=6
 for row in r['rows']:
  view=row['view'];ir=int(np.argmin(abs(np.array([.5,1,1.5])*7/.6766-row['radius_Mpc'])))
  full=d[f'{view}_{ir}_full'];nf=d[f'{view}_{ir}_nfull'];allowed=np.ones(n,bool) if view=='local' else d['original']
  projected=d[view+'_projected'];line=np.isfinite(projected)&(projected>=.5847730874520835)
  masks={'all':allowed,'true_wall':allowed&(t==1),'true_filament':allowed&(t==2),'predicted_wall':allowed&(p.argmax(1)==1)&(p[:,2]<.2),'line_predicted_wall':allowed&line&(p.argmax(1)==1)&(p[:,2]<.2),'line_true_wall':allowed&line&(t==1)}
  cohort=masks[row['cohort']];valid=cohort&d[f'{view}_{ir}_eligible']&(nf>=row['minimum'])&np.isfinite(full).all(1)
  assert cohort.sum()==row['cohort_rows'];assert np.array_equal(np.bincount(t[valid],minlength=4),row['truth_counts'])
  assert np.array_equal(np.bincount(full[valid].argmax(1),minlength=3),row['full_geometry_counts'])
  assert valid.sum()==row['wall_truth_fraction']['rows'];checks+=4
  if valid.any():
   for key,y in [('wall_truth_fraction',t==1),('filament_truth_fraction',t==2),('wall_probability',p[:,1]),('wall_calibration_residual',(t==1)-p[:,1]),('sheet_fraction',full.argmax(1)==1),('line_fraction',full.argmax(1)==0)]:
    np.testing.assert_allclose(row[key]['weighted_mean'],np.dot(w[valid],y[valid])/w[valid].sum(),rtol=1e-12,atol=1e-12);checks+=1
 # Verify all stored nonempty spectra sum to one and preserve count order.
 for view in ['local','original']:
  for ir in range(3):
   assert np.all(d[f'{view}_{ir}_nslab']<=d[f'{view}_{ir}_nfull']);checks+=1
   for key in ['full','slab','matched']:
    a=d[f'{view}_{ir}_{key}'];ok=np.isfinite(a).all(1)
    np.testing.assert_allclose(a[ok].sum(1),1,atol=1e-12);assert (a[ok]>=-1e-12).all();checks+=2
 report={'pass':True,'checks':checks,'rows':n,'geometry_sha256':digest(root/'geometry.npz'),'validator_sha256':digest(__file__)}
 (root/'VALIDATION.json').write_text(json.dumps(report,indent=2)+'\n')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axs=plt.subplots(1,2,figsize=(11,5),layout='constrained')
 rows=[a for a in r['rows'] if a['view']=='local' and a['minimum']==8 and abs(a['radius_Mpc']-7/.6766)<1e-8]
 get={a['cohort']:a for a in rows}
 for i,(key,label) in enumerate([('true_wall','True wall'),('true_filament','True filament'),('line_predicted_wall','String-like wall prediction')]):
  a=get[key]
  for j,(metric,color) in enumerate([('line_fraction','#a63603'),('sheet_fraction','#3182bd')]):
   v=a[metric]
   if v['rows']:
    axs[0].plot(i+(j-.5)*.18,v['weighted_mean'],'o',color=color,label=['3D line-dominant','3D sheet-dominant'][j] if i==0 else None)
    axs[0].vlines(i+(j-.5)*.18,*v['ci95'],color=color)
 axs[0].set_xticks(range(3),['True wall','True filament','String-like\nwall prediction']);axs[0].set(ylabel='Weighted fraction',ylim=(0,1),title='Full 3D tracer geometry');axs[0].legend(fontsize=8)
 for i,key in enumerate(['predicted_wall','line_predicted_wall']):
  a=get[key]
  for j,(metric,label,color) in enumerate([('wall_truth_fraction','Actual wall fraction','#3182bd'),('wall_probability','Mean P(wall)','#756bb1'),('filament_truth_fraction','Actual filament fraction','#d95f0e')]):
   v=a[metric]
   if v['rows']:
    axs[1].plot(i+(j-1)*.16,v['weighted_mean'],'o',color=color,label=label if i==0 else None);axs[1].vlines(i+(j-1)*.16,*v['ci95'],color=color)
 axs[1].set_xticks(range(2),['All low-Pfil\nwall predictions','String-like subset']);axs[1].set(ylim=(0,1),title='Tidal truth versus frozen predictions');axs[1].legend(fontsize=8)
 fig.suptitle('ph006 development mock · local 20 Mpc slabs · 0.20 ≤ z < 0.30\nR = 7 Mpc/h; ≥8 neighbours; conditional 95% block intervals',fontsize=12)
 fig.savefig(root/'mock_truth.png',dpi=170);fig.savefig(root/'mock_truth.pdf');print(json.dumps(report))
if __name__=='__main__':main()
