"""Compare producer branches individually with canonical parents and observed Loa."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mock_selection_test import digest,save

def run(root):
    repo=Path(__file__).resolve().parents[2]
    evidence=repo/'docs/evidence/p12a_alignment_execution_20260926'
    out=repo/'docs/figures/p12a_alignment_execution_20260926'
    loa=repo/'docs/evidence/p12a_joint_flags_20260926/HISTOGRAMS.npz'
    d=np.load(loa); a=np.load(evidence/'parent/v0.1.npz'); b=np.load(evidence/'parent/v1.npz')
    assert np.isclose(d['r_edges'][75],19.5)
    real=d['loa'][:,:,:75].sum(-1)
    z=(a['z_edges'][1:]+a['z_edges'][:-1])/2
    files={key:root/f'galaxy_cut_sky_{key}.npz' for key in ['N','S','N_Y3','S_Y3']}
    p={key:np.load(path) for key,path in files.items()}
    for t in p.values():assert np.array_equal(t['z_edges'],a['z_edges'])
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey='row',layout='constrained')
    for c,name in enumerate(['Galactic SGC','Galactic NGC']):
        series=[(a['z'][c],'canonical v0.1','#777777',':'),(b['z'][c],'canonical v1','#ac751c','--'),
                (p['N']['z'][c],'producer N','#2463a3','-'),(p['S']['z'][c],'producer S','#a3486c','-.')]
        for values,label,col,ls in series:
            ax[0,c].plot(z,values/d['area_deg2'][c]/.01,label=label,color=col,ls=ls)
        ax[0,c].plot(z,real[c]/d['area_deg2'][c]/.01,label='Loa observed',color='#222222',ls=(0,(5,2,1,2)))
        ax[0,c].set(title=name,yscale='log',ylabel='Counts / deg² / unit z')
        for branch,col,ls in [('N','#2463a3','-'),('S','#a3486c','-.')]:
            ax[1,c].plot(z,np.divide(p[branch+'_Y3']['z'][c],p[branch]['z'][c],out=np.full(40,np.nan),where=p[branch]['z'][c]>0),label=branch,color=col,ls=ls)
        ax[1,c].axhline(1,color='#222222',lw=1,ls=':')
        ax[1,c].set(xlabel='RSD redshift',ylabel='Y3 forFA / own raw counts')
    ax[0,0].legend(fontsize=9);ax[1,0].legend(title='Producer branch',fontsize=9)
    fig.suptitle('Exposed ph000 producer branches: common sky, 12 ≤ r < 19.5\nRaw/forFA and observed Loa are different stages; photometry parity unresolved',fontsize=12)
    fig.savefig(out/'producer_branches.png',dpi=160);fig.savefig(out/'producer_branches.pdf');plt.close(fig)
    save(out/'PRODUCER_PLOT_RECEIPT.json',dict(script_sha256=digest(__file__),
         inputs={str(q):digest(q) for q in [loa,evidence/'parent/v0.1.npz',evidence/'parent/v1.npz',*files.values()]},
         broad_counts={key:t['z'].reshape(2,4,10).sum(-1).tolist() for key,t in p.items()},
         loa_broad_counts=real.reshape(2,4,10).sum(-1).tolist(),
         outputs={q.name:digest(q) for q in [out/'producer_branches.png',out/'producer_branches.pdf']}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True)
    run(parser.parse_args().root)
