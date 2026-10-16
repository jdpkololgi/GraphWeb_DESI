"""Lightweight marginal plots from archived, audited 2-D histogram counts."""
from pathlib import Path
import csv
import hashlib
import json
import os
os.environ['MPLBACKEND'] = 'Agg'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / 'docs/evidence/loa_rest_colour_20261016'
DEST = ROOT / 'docs/figures/loa_colour_histograms_20261016'
EVID = ROOT / 'docs/evidence/loa_colour_histograms_20261016'
WINDOWS = [('z015020', '0.15 ≤ z < 0.20', 10.60),
           ('z020022', '0.20 ≤ z < 0.22', 10.69),
           ('z020030', '0.20 ≤ z < 0.30', 10.96)]
ENVS = ['Void', 'Wall', 'Filament', 'Cluster']
COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']
STYLES = ['-', '--', '-.', ':']

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    DEST.mkdir(exist_ok=True)
    EVID.mkdir(exist_ok=True)
    source = SRC / 'diagram_bins.npz'
    audit_path = SRC / 'diagram_validation.json'
    audit = json.loads(audit_path.read_text())
    bins = np.load(source, allow_pickle=False)
    edges = bins['colour_edges'][::2]
    widths = np.diff(edges)
    rows, checks = [], []
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    for method in ['hard', 'soft']:
        fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey='row')
        for col, (key, label, limit) in enumerate(WINDOWS):
            omitted = 0.
            for j, env in enumerate(ENVS):
                tag = f'observed_{method}_{key}_{j}'
                h = bins[tag]
                count = h.sum(axis=0).reshape(-1, 2).sum(axis=1)
                rec = audit[tag]
                assert np.isclose(count.sum(), rec['weight_displayed'], rtol=1e-12, atol=1e-7)
                assert np.all(count >= 0)
                density = count / count.sum() / widths
                assert np.isclose(np.sum(density * widths), 1, atol=1e-12)
                # Independent column-pair summation of the archived 2-D grid.
                reference = np.array([h[:, k:k+2].sum() for k in range(0, h.shape[1], 2)])
                assert np.allclose(count, reference, rtol=1e-12, atol=1e-7)
                omitted += rec['weight_outside']
                for row in range(2):
                    values = density if row == 0 else np.where(density > 0, density, np.nan)
                    axes[row, col].stairs(values, edges, color=COLORS[j],
                                          linestyle=STYLES[j], linewidth=1.8)
                checks.append(dict(method=method, window=key, environment=env,
                                   total=rec['weight_total'], displayed=float(count.sum()),
                                   omitted=rec['weight_outside'], integral=float(np.sum(density*widths))))
                for lo, hi, n, d in zip(edges[:-1], edges[1:], count, density):
                    rows.append(dict(method=method, window=key, environment=env,
                                     colour_lo=lo, colour_hi=hi, count_or_weight=n, density=d))
            axes[0, col].set_title(f'{label}\nlog₁₀ M★/M☉ ≥ {limit:.2f}', fontsize=12)
            axes[0, col].text(.03, .96, f'Outside saved bounds: {omitted:g}',
                              transform=axes[0, col].transAxes, va='top', fontsize=9)
            for row in range(2):
                ax = axes[row, col]
                ax.set_xlim(edges[0], edges[-1])
                ax.axvline(.75, color='#666666', linewidth=.8, linestyle='--', alpha=.6)
                ax.grid(axis='y', alpha=.18)
            axes[1, col].set_yscale('log')
            axes[1, col].set_ylim(1e-5, 10)
            axes[1, col].set_xlabel('Observed g − r [mag]')
        axes[0, 0].set_ylim(bottom=0)
        axes[0, 0].set_ylabel('Normalised density [mag⁻¹]')
        axes[1, 0].set_ylabel('Normalised density [mag⁻¹] · log scale')
        handles = [Line2D([], [], color=c, linestyle=s, lw=2, label=e)
                   for c, s, e in zip(COLORS, STYLES, ENVS)]
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .925), ncol=4, frameon=False)
        mode = 'Highest-probability environment' if method == 'hard' else 'Environment-probability weights'
        fig.suptitle(f'Observed colour distributions by inferred environment\n{mode} · original CIGALE mass selections', y=.99, fontsize=16)
        fig.subplots_adjust(left=.075, right=.985, top=.82, bottom=.14, hspace=.18, wspace=.12)
        fig.text(.075, .045,
                 'Each environment has unit area within the saved bounds; bin width = 0.0333 mag. Dashed vertical line: g − r = 0.75.\n'
                 'Lower panels repeat the same histograms on a log scale. Raw distributions: no mass/redshift matching or volume weighting.\n'
                 'Saved bounds: 10.3 ≤ log₁₀ M★/M☉ ≤ 12.5 and −0.5 ≤ g − r ≤ 2.5. Higher-redshift windows overlap.', fontsize=10)
        for ext in ['png', 'pdf']:
            fig.savefig(DEST / f'observed_colour_{method}.{ext}', dpi=170)
        plt.close(fig)
    with (EVID / 'histograms.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    manifest = dict(source_sha256=sha(source), audit_sha256=sha(audit_path),
                    script_sha256=sha(Path(__file__)), checks=checks,
                    normalization='Per-environment unit area within archived 2-D display bounds',
                    outputs={str(p.relative_to(ROOT)): sha(p) for p in sorted(DEST.iterdir())},
                    table_sha256=sha(EVID / 'histograms.csv'))
    (EVID / 'validation.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Validated {len(checks)} distributions; {len(rows)} CSV bins; output {DEST}')

if __name__ == '__main__':
    main()
