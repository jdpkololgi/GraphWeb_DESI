"""Recreate observed colour–mass diagrams from the three frozen selections."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
import socket
from pathlib import Path

BASE = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs')
OUT = BASE / 'loa_colour_mass_diagrams_20261016_v1'
CASES = [('z015020', 'loa_colour_lowz_20261016_v1', 'windows.json'),
         ('z020022', 'loa_colour_split_20261016_v2', 'windows.json'),
         ('z020030', 'loa_colour_split_20261016_v1', 'window_counts.json')]

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

def main():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Run inside an explicitly approved compute allocation')
    import fitsio
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    OUT.mkdir(exist_ok=False)
    # Shared bins and colour scale permit direct count comparisons.
    xe = np.linspace(10.5, 12.5, 161)
    ye = np.linspace(-.5, 2.5, 181)
    records, histograms = [], []
    for key, folder, windowfile in CASES:
        root = BASE / folder
        windows = json.loads((root / windowfile).read_text())
        w = windows.get('windows', windows)[key]
        path = root / ('LOA_COLOUR_MASS_COMPLETE_' + key + '.fits')
        a = fitsio.read(path, columns=['TARGETID', 'Z', 'MASS_CG_15', 'GR_OBS'])
        mass = np.log10(a['MASS_CG_15']); colour = a['GR_OBS']
        zlo = w.get('zlo', .2)
        assert len(a) == w['rows'] and len(np.unique(a['TARGETID'])) == len(a)
        assert ((a['Z'] >= zlo) & (a['Z'] < w['zhi']) & (mass >= w['logmass'])).all()
        valid = np.isfinite(mass) & np.isfinite(colour)
        red = valid & (colour > .75) & (colour < 1.5)
        blue = valid & (colour < .75)
        assert red.sum() == w['counts']['red'] and blue.sum() == w['counts']['blue']
        h, _, _ = np.histogram2d(mass[valid], colour[valid], bins=(xe, ye))
        displayed = valid & (mass >= xe[0]) & (mass <= xe[-1]) & (colour >= ye[0]) & (colour <= ye[-1])
        assert int(h.sum()) == displayed.sum()
        # Independent flat-bin count cross-check, including the final edges.
        xi = np.minimum(np.searchsorted(xe, mass[displayed], side='right')-1, len(xe)-2)
        yi = np.minimum(np.searchsorted(ye, colour[displayed], side='right')-1, len(ye)-2)
        check = np.bincount(xi*(len(ye)-1)+yi, minlength=h.size).reshape(h.shape)
        assert np.array_equal(h, check)
        histograms.append(h)
        records.append(dict(key=key, zlo=zlo, zhi=w['zhi'], logmass=w['logmass'],
            rows=len(a), red=int(red.sum()), blue=int(blue.sum()), excluded=int((~(red|blue)).sum()),
            nonfinite=int((~valid).sum()), outside_axes=int((valid & ~displayed).sum()),
            displayed=int(displayed.sum()), source=str(path), source_sha256=sha(path),
            windows_sha256=sha(root/windowfile)))
    plt.rcParams.update({'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False})
    norm = LogNorm(vmin=1, vmax=max(h.max() for h in histograms))
    def draw(ax, rec, hist):
        mesh=ax.pcolormesh(xe, ye, np.ma.masked_equal(hist.T,0), norm=norm, cmap='magma', rasterized=True)
        ax.axvspan(xe[0],rec['logmass'],color='#e6e6e6',zorder=0)
        ax.axvline(rec['logmass'],color='#222222',ls='--',lw=1.3)
        ax.axhline(.75,color='#168ac1',ls='--',lw=1.2)
        ax.axhline(1.5,color='#b43d4b',ls=':',lw=1.2)
        ax.set(xlim=(xe[0],xe[-1]),ylim=(ye[0],ye[-1]),xlabel=r'$\log_{10}(M_\star/M_\odot)$')
        ax.set_title(f"{rec['zlo']:.2f} ≤ z < {rec['zhi']:.2f}  |  log M★ ≥ {rec['logmass']:.2f}\nRed {rec['red']:,} · blue {rec['blue']:,}",fontsize=11)
        ax.text(.03,.97,f"N = {rec['rows']:,}\nOutside colour cuts: {rec['excluded']:,}\nOutside axes: {rec['outside_axes']:,}; nonfinite: {rec['nonfinite']:,}",transform=ax.transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',alpha=.88,edgecolor='none'))
        return mesh
    fig, axs = plt.subplots(1,3,figsize=(17,6),sharex=True,sharey=True)
    fig.subplots_adjust(left=.06,right=.90,bottom=.20,top=.80,wspace=.13)
    for ax,r,h in zip(axs,records,histograms):mesh=draw(ax,r,h)
    axs[0].set_ylabel('Observed g − r [mag; Milky-Way corrected]')
    fig.colorbar(mesh,cax=fig.add_axes([.92,.20,.012,.60]),label='Galaxies per 0.0125 dex × 0.0167 mag bin')
    fig.suptitle('Loa colour–mass distributions · fixed empirical completeness selections',y=.96,fontsize=16)
    fig.text(.06,.07,'Red: 0.75 < g−r < 1.5; blue: g−r < 0.75. All selected galaxies shown, including those outside these cuts.\nNo K correction. Shared logarithmic count scale; raw populations, without environment weighting or mass/redshift matching.',fontsize=10)
    for ext in ['png','pdf']:fig.savefig(OUT/f'colour_mass_comparison.{ext}',dpi=180)
    plt.close(fig)
    for r,h in zip(records,histograms):
        fig,ax=plt.subplots(figsize=(7,6));mesh=draw(ax,r,h)
        ax.set_ylabel('Observed g − r [mag; Milky-Way corrected]')
        fig.colorbar(mesh,ax=ax,label='Galaxies per bin (shared scale)')
        fig.text(.08,.02,'Blue: g−r < 0.75; red: 0.75 < g−r < 1.5. No K correction.',fontsize=9)
        fig.tight_layout(rect=[0,.05,1,1])
        for ext in ['png','pdf']:fig.savefig(OUT/f"colour_mass_{r['key']}.{ext}",dpi=180)
        plt.close(fig)
    np.savez_compressed(OUT/'histograms.npz',mass_edges=xe,colour_edges=ye,**{r['key']:h for r,h in zip(records,histograms)})
    manifest=dict(job=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),code_sha256=sha(Path(__file__)),selections=records,
        validation='Exact saved counts, unique IDs, redshift and mass limits, colour membership, histogram sum and independent bin reconstruction passed',
        outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file()})
    (OUT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2),flush=True)

if __name__ == '__main__':
    main()
