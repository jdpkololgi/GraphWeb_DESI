"""Release and angular-boundary sensitivity figures from completed receipts."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mock_selection_test import digest,save

def ratio(a,b):return np.divide(a,b,out=np.full(a.shape,np.nan),where=b>0)
def run(dr1,boundary,out):
    out.mkdir(parents=True,exist_ok=True)
    r=np.load(dr1/'HISTOGRAMS.npz');b=np.load(boundary/'HISTOGRAMS.npz');assert np.array_equal(r['z_edges'],b['z_edges'])
    z=(r['z_edges'][1:]+r['z_edges'][:-1])/2
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey='row',layout='constrained')
    for c,name in enumerate(['SGC','NGC']):
        for key,label,col,ls in [('original','Original mask: includes non-DR1 sky','#777777',':'),('random256','Both randoms, NSIDE256','#2463a3','--'),('random512','Both randoms, NSIDE512','#a3486c','-.'),('random512_interior','Both randoms, interior','#ac751c','-')]:
            ax[0,c].plot(z,ratio(r['loa_'+key][c],r['iron_'+key][c]),label=label,color=col,ls=ls)
        ax[0,c].set(title=name,ylabel='Loa / DR1 observed counts')
        for key,label,col,ls in [('original','Original support','#777777',':'),('erode1','Remove 1 pixel ring','#2463a3','--'),('erode2','Remove 2 pixel rings','#ac751c','-')]:
            ax[1,c].plot(z,ratio(b['loa_'+key][c].sum(0),b['mock_'+key][c].sum(0)),label=label,color=col,ls=ls)
        ax[1,c].set(xlabel='Redshift',ylabel='Loa / original observed mock')
        for a in ax[:,c]:a.axhline(1,color='#222222',lw=1,ls=':')
    ax[0,0].legend(fontsize=9);ax[1,0].legend(fontsize=9)
    fig.suptitle('Observed-stage count controls: release and footprint sensitivity\nFinite-random support and pixel erosion are diagnostics, not exact survey masks',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(out/f'release_boundary_controls.{ext}',dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),sharey=True,layout='constrained')
    for c,name in enumerate(['SGC','NGC']):
        for nt,label,col,ls in [(0,'NTILE ≤ 1','#2463a3','-'),(1,'NTILE = 2','#ac751c','--'),(2,'NTILE ≥ 3','#a3486c','-.')]:
            ax[c].plot(z,ratio(b['loa_erode2'][c,nt],b['mock_erode2'][c,nt]),label=label,color=col,ls=ls)
        ax[c].axhline(1,color='#222222',ls=':',lw=1);ax[c].set(title=name,xlabel='Redshift',ylabel='Loa / original observed mock')
    ax[0].legend(fontsize=9);fig.suptitle('Two-ring interior: tile-count strata\nSame NTILE does not guarantee equal assignment completeness',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(out/f'tile_strata_controls.{ext}',dpi=160)
    plt.close(fig)
    save(out/'PLOT_RECEIPT.json',dict(script_sha256=digest(__file__),inputs={str(p):digest(p) for p in [dr1/'HISTOGRAMS.npz',boundary/'HISTOGRAMS.npz']},outputs={p.name:digest(p) for p in out.glob('*.png')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dr1',type=Path,required=True);p.add_argument('--boundary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.dr1,a.boundary,a.out)
