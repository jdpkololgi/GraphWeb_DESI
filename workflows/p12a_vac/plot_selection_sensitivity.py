"""Summarise paired perturbation results; selected-core diagnostics only."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_selection_sensitivity_20260925_v1')
OUT=Path(__file__).resolve().parents[2]/'docs/figures/p12a_selection_20260925'
def main():
 result=json.loads((ROOT/'RESULTS.json').read_text());cores=[r['core'] for r in result['runs'] if r['arm']=='baseline'];arms=['baseline','half_expected','thin_wrong_expected','thin_matched_expected'];labels=['Baseline','Same counts\n½ expected','½ counts\noriginal expected','½ counts\n½ expected'];summary={};basep=None;baseq=None;fig,ax=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
 for i,arm in enumerate(arms):
  d=[np.load(ROOT/f'core{core}_{arm}.npz') for core in cores];p=np.concatenate([x['p_web'] for x in d]);q=np.concatenate([x['quantiles'] for x in d]);truth=np.concatenate([x['truth'] for x in d]);tc=(truth>.2).sum(-1)
  for core,x in zip(cores,d):assert np.array_equal(x['parent_id'],np.load(ROOT/f'core{core}_baseline.npz')['parent_id'])
  if arm=='baseline':basep=p;baseq=q
  c68=((truth>=q[:,1])&(truth<=q[:,3])).mean(0);c90=((truth>=q[:,0])&(truth<=q[:,4])).mean(0);r=dict(rows=len(p),mean_probability=p.mean(0).tolist(),coverage68=c68.tolist(),coverage90=c90.tolist(),mean_probability_shift=(p-basep).mean(0).tolist(),brier=float(((p-np.eye(4)[tc])**2).sum(1).mean()),median_shift_over_baseline_width=np.median((q[:,2]-baseq[:,2])/np.maximum(baseq[:,3]-baseq[:,1],1e-6),axis=0).tolist());summary[arm]=r
  for j,col in enumerate(['#4477aa','#ccbb44','#228833','#aa3377']):ax[0].bar(i+(j-1.5)*.18,r['mean_probability_shift'][j],width=.18,color=col,label=['Void','Sheet','Filament','Knot'][j] if i==0 else None)
  for j in range(3):ax[1].scatter(i+(j-1)*.1,c68[j],marker=['o','s','^'][j],color=['#0072b2','#d55e00','#009e73'][j],label=f'λ{j+1}' if i==0 else None)
  ax[2].bar(i,r['brier'],color='#4477aa')
 for a in ax:a.set_xticks(range(4),labels,fontsize=9)
 ax[0].axhline(0,color='grey',lw=1);ax[0].set_ylabel('Change in mean class probability');ax[0].legend(fontsize=8,ncol=2);ax[1].axhline(.68,color='black',ls=':',label='Nominal 68%');ax[1].set_ylabel('Empirical central 68% coverage');ax[1].legend(fontsize=8);ax[2].set_ylabel('Mean multiclass Brier score (lower is better)')
 fig.suptitle('Density perturbations: matching the mean does not restore the original information',fontsize=15);fig.supxlabel(f'{len(p)} identical retained query galaxies, eight preselected ph006 cores, 512 draws each; paired RNG.\nFixed angular response. Local development diagnostic, not survey-wide coverage or an approved DESI correction.',fontsize=9)
 for ext in ['png','pdf']:fig.savefig(OUT/f'04_selection_sensitivity.{ext}',dpi=170,bbox_inches='tight')
 (OUT/'SENSITIVITY.json').write_text(json.dumps(dict(summary=summary,source=str(ROOT/'RESULTS.json'),baseline_replay_pass=result['baseline_replay_pass']),indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
