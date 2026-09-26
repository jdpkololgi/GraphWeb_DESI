"""Render audited count/version screens from compact census products."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mock_selection_test import digest,save

def run(parent,uchuu,out,internal=None):
    out.mkdir(parents=True,exist_ok=True)
    a=np.load(parent/'v0.1.npz');b=np.load(parent/'v1.npz')
    loa=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_joint_flags_20260926/HISTOGRAMS.npz'
    d=np.load(loa);z=(a['z_edges'][1:]+a['z_edges'][:-1])/2
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(11,7),sharex=True,layout='constrained')
    for c,name in enumerate(['SGC','NGC']):
        for x,label,col,ls in [(a['z'][c],'v0.1 intrinsic','#2463a3','--'),(b['z'][c],'v1 intrinsic','#b77419','-')]:
            ax[0,c].plot(z,x/d['area_deg2'][c]/.01,label=label,color=col,ls=ls)
        # Same numeric r cut, current baseline quality. No completeness correction.
        real=d['loa'][c,:,:75].sum(-1)
        ax[0,c].plot(z,real/d['area_deg2'][c]/.01,color='#333333',label='Loa observed, 12≤r<19.5',ls=':')
        ax[0,c].set(title=name,yscale='log',ylabel='Counts / deg² / unit z')
        ax[1,c].plot(z,b['z'][c]/a['z'][c],color='#2463a3')
        ax[1,c].axhline(1,color='#333333',ls='--',lw=1)
        ax[1,c].set(xlabel='Redshift',ylabel='Intrinsic v1 / v0.1',ylim=(.7,1.08*np.max(b['z']/a['z'])))
    ax[0,0].legend(fontsize=9)
    fig.suptitle('Abacus ph000: same sky and numerical magnitude cut\nIntrinsic vs observed curves are different selection stages',fontsize=13)
    fig.savefig(out/'abacus_versions.png',dpi=160);fig.savefig(out/'abacus_versions.pdf');plt.close(fig)
    fig,ax=plt.subplots(2,3,figsize=(12,7),sharex=True,sharey=True,layout='constrained')
    for c,name in enumerate(['SGC','NGC']):
        for j,(h,title) in enumerate([(a['zr'][c],'v0.1 intrinsic'),(b['zr'][c],'v1 intrinsic'),(d['loa'][c],'Loa observed')]):
            density=h[:,:75]/d['area_deg2'][c]/(.01*.1)
            mesh=ax[c,j].pcolormesh(a['z_edges'],a['r_edges'][:76],np.ma.masked_equal(density.T,0),norm=LogNorm(.1,5000),cmap='cividis',rasterized=True)
            ax[c,j].set(title=f'{name}: {title}',xlabel='Redshift',ylabel='Apparent r',ylim=(19.5,15))
    fig.colorbar(mesh,ax=ax,label='Counts / deg² / dz / dr',shrink=.85)
    fig.suptitle('Magnitude–redshift census; passband equivalence remains unproven')
    fig.savefig(out/'magnitude_redshift.png',dpi=160);fig.savefig(out/'magnitude_redshift.pdf');plt.close(fig)
    u=np.load(uchuu/'HISTOGRAMS.npz');z=(u['z_edges'][1:]+u['z_edges'][:-1])/2
    ratios=u['counts']/95.7/u['data_per_deg2'][None,...]
    fig,ax=plt.subplots(1,2,figsize=(11,4),sharey=True,layout='constrained')
    for c,name in enumerate(['S','N']):
        mean=ratios[:,c].mean(0);std=ratios[:,c].std(0,ddof=1)
        ax[c].plot(z,mean,color='#2463a3',label='102-realization mean')
        ax[c].fill_between(z,mean-std,mean+std,color='#2463a3',alpha=.2,label='±1 realization SD')
        ax[c].plot(z,u['fuji31_per_deg2'][c]/u['data_per_deg2'][c],color='#b77419',ls=':',label='Fuji3.1 / archived data')
        ax[c].axhline(1,color='#333333',ls='--',lw=1)
        ax[c].set(title=name,xlabel='Redshift',ylabel='Area-normalized count ratio',ylim=(0,1.7))
    ax[0].legend(fontsize=9)
    fig.suptitle('Uchuu public SV3 / archived SV3 data\nPublished Nbin treatment inherited; this is not a Loa qualification',fontsize=12)
    fig.savefig(out/'uchuu_sv3.png',dpi=160);fig.savefig(out/'uchuu_sv3.pdf');plt.close(fig)
    inputs=[parent/'v0.1.npz',parent/'v1.npz',uchuu/'HISTOGRAMS.npz',loa]
    assert np.isclose(d['r_edges'][75],19.5)
    if internal is not None:
        inputs.append(internal/'HISTOGRAMS.npz')
        t=np.load(internal/'HISTOGRAMS.npz');zz=(t['z_edges'][1:]+t['z_edges'][:-1])/2
        fig,ax=plt.subplots(1,2,figsize=(11,4),sharey=True,layout='constrained')
        for c,name in enumerate(['SGC','NGC']):
            for key,label,col,ls in [('glam_bgs_v2','GLAM BGS v2, mock10','#2463a3','--'),('holi_bgs_v2','Holi BGS v2, mock0','#b77419','-')]:
                ax[c].plot(zz,t[key][c]/t['loa'][c],color=col,ls=ls,label=label)
            ax[c].axhline(1,color='#333333',ls=':',lw=1)
            ax[c].axvline(.5,color='#777777',ls=':',lw=1)
            ax[c].set(title=name,xlabel='Redshift',ylabel='Observed-stage mock / Loa counts',ylim=(-.05,1.3))
        ax[0].legend(fontsize=9)
        fig.suptitle('Internal Loa mock screen on shared occupied sky\nOne realization each; no DELTACHI2/SPECTYPE counterparts',fontsize=12)
        fig.savefig(out/'internal_loa.png',dpi=160);fig.savefig(out/'internal_loa.pdf');plt.close(fig)
    save(out/'PLOT_RECEIPT.json',dict(script_sha256=digest(__file__),inputs={str(p):digest(p) for p in inputs},
        outputs={p.name:digest(p) for p in out.glob('*.png')},
        parent_broad_ratio=(b['z'].reshape(2,4,10).sum(-1)/a['z'].reshape(2,4,10).sum(-1)).tolist(),
        parent_central_fraction={v:(x['central'].reshape(2,4,10).sum(-1)/x['z'].reshape(2,4,10).sum(-1)).tolist() for v,x in [('v0.1',a),('v1',b)]}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,required=True);p.add_argument('--uchuu',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--internal',type=Path)
    a=p.parse_args();run(a.parent,a.uchuu,a.out,a.internal)
