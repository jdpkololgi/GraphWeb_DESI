"""Internal Uchuu count ratios with distinct data-selection/weight controls."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mock_selection_test import digest,save

def run(parent,parity,out):
    a=np.load(parent/'HISTOGRAMS.npz');b=np.load(parity/'HISTOGRAMS.npz')
    assert np.array_equal(a['uchuu_altmtl'],b['mock'])
    z=(a['z_edges'][1:]+a['z_edges'][:-1])/2
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(11,4.5),sharey=True,layout='constrained')
    ratios={}
    for c,name in enumerate(['SGC','NGC']):
        series=[(b['mock'][c]/b['current_loa'][c],'Raw mock / current Loa full','#2463a3','-'),
                (a['uchuu_altmtl'][c]/a['paired_loa'][c],'Raw mock / archived Loa clustering','#b77419','--'),
                (a['uchuu_altmtl_weighted'][c]/a['paired_loa_weighted'][c],'Stored-weight mock / archived data','#555555',':')]
        for y,label,color,style in series:ax[c].plot(z,y,label=label,color=color,ls=style)
        ax[c].axhline(1,color='#333333',lw=.8)
        ax[c].set(title=name,xlabel='Redshift',ylabel='Mock / data count ratio')
        ratios[name]={label:y.tolist() for y,label,_,_ in series}
    ax[0].legend(fontsize=8,loc='lower left')
    fig.suptitle('Internal Uchuu Y3-v2.0 BGS: selection and weighting controls\nOne realization; exact observation and photometry parity remains open',fontsize=12)
    out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/'uchuu_dr2_loa.png',dpi=160);fig.savefig(out/'uchuu_dr2_loa.pdf');plt.close(fig)
    save(out/'UCHUU_DR2_PLOT_RECEIPT.json',dict(script_sha256=digest(__file__),
        inputs={str(p):digest(p) for p in [parent/'HISTOGRAMS.npz',parity/'HISTOGRAMS.npz']},
        outputs={name:digest(out/name) for name in ['uchuu_dr2_loa.png','uchuu_dr2_loa.pdf']},ratios=ratios))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,required=True);p.add_argument('--parity',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.parent,a.parity,a.out)
