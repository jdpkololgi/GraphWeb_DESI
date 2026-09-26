"""Standalone scientific diagnostic figure for the isolated ph006 selection test."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    s=json.loads((a.root/'SWEEP.json').read_text());prep=json.loads((a.root/'PREPARED.json').read_text())
    base=next(x for x in s['trials'] if x['delta_m']==0 and x['extra_q']==0);best=s['best'];z=np.linspace(.155,.545,40);target=np.array(s['target_counts'])
    review=[]
    for row in s['trials']:
        ratio=np.array(row['fixed_retention_prediction'])/target
        review.append(dict(delta_m=row['delta_m'],extra_q=row['extra_q'],sgc_fine_log_mse=float(np.mean(np.log(ratio[0])**2)),fine_ratio_min=ratio.min(1).tolist(),fine_ratio_max=ratio.max(1).tolist(),max_fractional_error_by_cap=np.max(abs(ratio-1),axis=1).tolist()))
    fine=dict(best_existing_grid_by_fine_sgc=min(review,key=lambda x:x['sgc_fine_log_mse']),any_existing_trial_with_all_bins_within_10_percent=any(max(x['max_fractional_error_by_cap'])<=.1 for x in review),fine_bin_width=.01,posthoc_diagnostic=True,trials=review)
    (a.root/'FINE_BIN_REVIEW.json').write_text(json.dumps(fine,indent=2)+'\n')
    fig,ax=plt.subplots(2,2,figsize=(10,7),sharex=True)
    for cap,label in enumerate(['SGC: tuning region','NGC: exposed check region']):
        for row,name,style in [(base,'Baseline','--'),(best,'Best fixed-retention screen','-')]:
            ratio=np.array(row['fixed_retention_prediction'])[cap]/target[cap]
            ax[0,cap].plot(z,ratio,style,label=name,lw=1.6)
        q40=np.array(prep['stages']['quality40_r19p5'])[cap]/target[cap]
        ax[0,cap].plot(z,q40,':',color='grey',label='Loa quality40 / quality25')
        ax[0,cap].axhline(1,color='black',lw=.7);ax[0,cap].set_title(label);ax[0,cap].set_ylim(0,1.6)
        ax[0,cap].set_ylabel('Counts / Loa quality25')
        for key,name,style in [('baseline','Baseline','--'),('screening_best','Screening candidate','-')]:
            rows=[x for x in s['population'][key] if x['cap']==cap]
            ax[1,cap].plot([np.mean(x['z']) for x in rows],[x['mock_colour_median']-x['real_colour_median'] for x in rows],style+'o',label=name)
        ax[1,cap].axhline(0,color='black',lw=.7);ax[1,cap].set_ylabel('Mock − Loa median g−r [mag]');ax[1,cap].set_xlabel('Redshift');ax[1,cap].set_ylim(-.3,.1)
        for i in range(2):ax[i,cap].grid(alpha=.2)
    ax[0,0].legend(fontsize=8);ax[1,0].legend(fontsize=8)
    fig.suptitle('ph006 development test: count agreement is not population validation')
    fig.text(.5,.015,'Counts use fixed baseline retention, not fresh fibre assignment. Colour compares selected raw mocks with successful Loa; passbands remain unresolved.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.04,1,.95]);fig.savefig(a.root/'selection_screen.png',dpi=160);fig.savefig(a.root/'selection_screen.pdf');plt.close(fig)
if __name__=='__main__':main()
