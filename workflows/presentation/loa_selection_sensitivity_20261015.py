"""Fixed-prediction selection sensitivity; explicitly the MLP, not GAT.

Three declared redshift windows and three completeness choices. Never choose
a preferred model/cut from the resulting observational trend strength.
"""
import json
from pathlib import Path
import numpy as np
from loa_tng_abacus_overlay_20261015 import (
    HOME, OUT as INPUT, KEYS, LABELS, PROPLABELS, property_data,
    memberships, raw_statistics, sha, dump)
from legacy_selection_tests_20261014 import completeness, apply_limit
from desi_density_properties_20261014 import estimate

OUT = HOME/'loa_selection_sensitivity_20261015'
WINDOWS = [('full', .2, .3), ('lower', .2, .25), ('upper', .25, .3)]


def run():
    OUT.mkdir(exist_ok=True)
    source = INPUT/'loa_predictions.npz'
    assert sha(source) == json.loads((INPUT/'inference_receipt.json').read_text())['prediction_sha256']
    a = np.load(source)
    photpath = HOME/'legacy_selection_tests_20261014/p12a_targeting.npz'
    phot = np.load(photpath)
    assert np.array_equal(a['TARGETID'], phot['TARGETID'])
    good, yy = property_data(a)
    bright = phot['r'] < phot['rlimit']
    parent = a['CIGALE_MATCHED'] & (a['MASS_CG_15'] > 0) & np.isfinite(yy['mass']) & bright
    edges = np.linspace(.2, .3, 11)
    curve = completeness(a['z'][parent], yy['mass'][parent], phot['r'][parent], phot['rlimit'][parent], edges)
    cuts = {'uncut':np.ones(len(good), bool),
            'complete90':apply_limit(a['z'], yy['mass'], curve),
            'complete90_plus0p2':apply_limit(a['z'], yy['mass'], curve, .2)}
    # Freeze group definitions on the complete valid parent; cuts do not relabel.
    base = good & bright
    ix = np.flatnonzero(base)
    zb = np.floor((a['z'][ix]-.2+1e-10)/.01).astype(int)
    rank = {k:memberships(a[k][ix].sum(1), zb) for k in KEYS}
    classes = {k:np.eye(4)[(a[k][ix] >= .2).sum(1)] for k in KEYS}
    result = {'curve':curve, 'windows':WINDOWS, 'scenarios':{}, 'paired_cuts':{},
              'groups':'Fixed per-model trace quartiles in dz=.01 on valid uncut bright parent; physical classes count mean eigenvalues >=.2.',
              'models':'TNG geometry-only MLPs, three-seed mean; Abacus saved P12-A. NOT the old GAT.',
              'limits':'Exploratory selection sensitivity on available .2-.3 predictions; no low-z inference or recovery claim from selecting the strongest panel. Completeness is a flux-limit diagnostic, not correction for all survey selection.',
              'inputs':{str(p):sha(p) for p in [source, photpath]}, 'verification':{}}
    checks = 0
    for window, lo, hi in WINDOWS:
        wz = (a['z'][ix] >= lo) & (a['z'][ix] < hi)
        for cut, keep in cuts.items():
            take = wz & keep[ix]
            rows = ix[take]
            _, region = np.unique(a['region'][rows], return_inverse=True)
            nr = region.max()+1
            draws = np.random.default_rng(20261015).multinomial(nr, np.ones(nr)/nr, size=128)
            weights = {k:v[take] for k,v in rank.items()}
            cl = {k:v[take] for k,v in classes.items()}
            mb = np.floor((yy['mass'][rows]-9+1e-10)/.1).astype(int)
            suites = {'rank_raw':{}, 'rank_controlled':{}, 'class_raw':{}}
            for prop, values in yy.items():
                y = values[rows].astype(float)
                suites['rank_raw'][prop] = raw_statistics(y, weights, region, draws)
                suites['class_raw'][prop] = raw_statistics(y, cl, region, draws)
                cells = zb[take] if prop == 'mass' else zb[take]+11*mb
                try:
                    r, common = estimate(y, weights, cells, region, np.zeros(len(y)), minimum=20, nboot=128)
                    r.pop('paired_density_minus_soft_contrast', None)
                    suites['rank_controlled'][prop] = r
                    for k,w in weights.items():
                        means = np.zeros(4)
                        for cell, ref in zip(r['common_cell_ids'], r['reference_weights']):
                            cc = cells == cell
                            means += ref*(w[cc].T@y[cc])/w[cc].sum(0)
                        assert np.allclose(means,r['methods'][k]['mean'],atol=1e-10,rtol=0)
                        checks += 1
                except ValueError as e:
                    suites['rank_controlled'][prop] = {'unavailable':str(e)}
                for k,w in weights.items():
                    expected = (w.T@y)/w.sum(0)
                    assert np.allclose(expected,suites['rank_raw'][prop]['methods'][k]['mean'],atol=1e-10,rtol=0)
                    checks += 1
            result['scenarios'][window+'_'+cut] = {'rows':len(rows), 'z':[lo,hi], 'cut':cut, 'suites':suites}
            print(window,cut,len(rows),{k:round(v['contrast']*100,2) for k,v in suites['rank_controlled']['quenched'].get('methods',{}).items()},flush=True)
        # A matched parent/cut reference mixture distinguishes changing mass mix
        # from changing environment association. Same rank definitions throughout.
        rows = ix[wz]
        cells = zb[wz]+11*np.floor((yy['mass'][rows]-9+1e-10)/.1).astype(int)
        _, region = np.unique(a['region'][rows], return_inverse=True)
        for cut in ['complete90','complete90_plus0p2']:
            weights = {k+'_parent':v[wz] for k,v in rank.items()}
            weights.update({k+'_cut':v[wz]*cuts[cut][rows,None] for k,v in rank.items()})
            try:
                r,_ = estimate(yy['quenched'][rows],weights,cells,region,np.zeros(len(rows)),minimum=20,nboot=128)
                r.pop('paired_density_minus_soft_contrast',None)
                for k,w in weights.items():
                    means = np.zeros(4)
                    for cell,ref in zip(r['common_cell_ids'],r['reference_weights']):
                        cc = cells == cell
                        means += ref*(w[cc].T@yy['quenched'][rows][cc])/w[cc].sum(0)
                    assert np.allclose(means,r['methods'][k]['mean'],atol=1e-10,rtol=0)
                    checks += 1
            except ValueError as e:
                r = {'unavailable':str(e)}
            result['paired_cuts'][window+'_'+cut] = r
    result['verification'] = {'independent_mean_checks':checks, 'targetid_parity':True,
        'script_sha256':sha(__file__), 'nested_cuts':bool(np.all(~cuts['complete90_plus0p2']|cuts['complete90']))}
    assert result['verification']['nested_cuts']
    dump(OUT/'RESULTS.json', result)
    render(result)


def render(r):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from shared.plot_style import apply_style, ACCENT_COLORS, TEXT_COLOR
    apply_style()
    colors = [TEXT_COLOR, ACCENT_COLORS['blue'], ACCENT_COLORS['blue'], ACCENT_COLORS['magenta'], ACCENT_COLORS['magenta']]
    with PdfPages(OUT/'selection_overlays.pdf') as pdf:
        for suite in ['rank_raw','rank_controlled','class_raw']:
            for prop in PROPLABELS:
                f,axs = plt.subplots(3,3,figsize=(16,13),sharex=True,sharey=True)
                for row,(window,lo,hi) in enumerate(WINDOWS):
                    for col,cut in enumerate(['uncut','complete90','complete90_plus0p2']):
                        ax=axs[row,col];s=r['scenarios'][window+'_'+cut];p=s['suites'][suite][prop]
                        for j,k in enumerate(KEYS):
                            if 'methods' not in p:continue
                            d=p['methods'][k];mu=np.asarray(d['mean'],float);ci=np.asarray(d['interval'],float)
                            x=np.arange(4)+(j-2)*.025
                            ax.plot(x,mu,marker='o',ms=3,color=colors[j],ls='--' if k.endswith('dense') else '-',label=LABELS[k])
                            ax.vlines(x,ci[0],ci[1],color=colors[j],alpha=.6)
                        ax.set_title(f'{lo:.2f}–{hi:.2f} · {cut}\nN={p.get("rows",s["rows"]):,}',fontsize=11)
                        ax.set_xticks(range(4),['Void','Wall','Filament','Cluster'] if suite=='class_raw' else ['Q1','Q2','Q3','Q4'])
                        ax.grid(alpha=.15)
                        if col==0:ax.set_ylabel(PROPLABELS[prop])
                        if 'methods' not in p:ax.text(.5,.5,'Insufficient common support',ha='center',transform=ax.transAxes)
                f.suptitle(f'Selection sensitivity: {prop} · {suite}\nTNG local-feature MLPs and Abacus P12-A; same Loa galaxies per panel',fontsize=17)
                handles,labels=axs[0,0].get_legend_handles_labels()
                f.legend(handles,labels,loc='lower center',ncol=3,bbox_to_anchor=(.5,.035))
                f.text(.5,.012,'Frozen rank groups; 128 sky-block draws, 16–84%. Cuts change population; controlled panels use mass/redshift cells.',ha='center',fontsize=10)
                f.subplots_adjust(top=.90,bottom=.12,hspace=.32,wspace=.12)
                pdf.savefig(f)
                if prop in ['quenched','ssfr']:
                    f.savefig(OUT/f'{suite}_{prop}.png',dpi=130)
                plt.close(f)


if __name__=='__main__':run()
