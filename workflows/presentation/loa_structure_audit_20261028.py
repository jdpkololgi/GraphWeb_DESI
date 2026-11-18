"""Bounded saved-field diagnostics. No inference, truth access, or release overwrite."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket

import h5py
import numpy as np
from scipy import fft, ndimage

SOURCE = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1/production_v1/NGC_seed17/unet_x_17_FIELD.h5')
OPERATOR = Path('/global/u2/d/dkololgi/TNG/Illustris/workflows/sbi/e2e_footprint_tidal_operator.py')
EXPECTED = 'b5e6f7173ce69d0146e7441004d5553c3822ebd042b44a5fdb26f2104db30b7e'
CELL = 6.766
CENTER = np.array([-554.6517793611433, -236.19830282255862, 541.6343589342335])


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main(out):
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Approved compute allocation required')
    out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('tidal', OPERATOR)
    op = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(op)
    assert sha(SOURCE) == EXPECTED
    checks = {}
    # Independent direct summation, translation/padding and contact-term checks.
    rng = np.random.default_rng(98)
    x = rng.normal(size=(5, 6, 7))
    target = (2, 3, 2)
    offset = np.array([1, 2, 3])
    a = op.convolve(x, op.kernel(x.shape, target, offset), target, workers=4)
    b = np.zeros_like(a)
    for index in np.ndindex(target):
        for source_index in np.ndindex(x.shape):
            b[index] += x[source_index] * op.cube_tensor(*(np.array(index)+offset-source_index))
    checks['direct_sum_maxabs'] = float(abs(a-b).max())
    padded = op.convolve(x, op.kernel(x.shape, target, offset, extra_padding=7), target, workers=4)
    checks['padding_maxabs'] = float(abs(a-padded).max())
    checks['self_term_maxabs'] = float(abs(op.cube_tensor(0.,0.,0.)-np.array([1/3,0,0,1/3,0,1/3])).max())
    assert max(checks.values()) < 1e-10
    results = []
    layers = {}
    with h5py.File(SOURCE) as f:
        origin = np.asarray(f.attrs['origin_mpc_h']) + CELL/2
        center = np.rint((CENTER-origin)/CELL).astype(int)
        assert np.max(abs(origin+center*CELL-CENTER)) < 1e-8
        target_shape = (32, 32, 32)
        target_lo = center-16
        tcrop = tuple(slice(int(v), int(v+32)) for v in target_lo)
        support_core = f['support'][tcrop]
        coverage_core = f['coverage'][tcrop]
        core_valid = (support_core >= .5) & (coverage_core > 0)
        assert core_valid.any(), 'No supported targets in selected core'
        core_rho = f['rho'][(slice(None), *tcrop)].astype('f8')
        # Compare fixed central target with increasing source context. These are
        # sensitivity probes, not convergence to an unknown cosmological exterior.
        for size in (96, 144, 192):
            lo = center-size//2
            assert np.all(lo >= 0) and np.all(lo+size <= f['coverage'].shape)
            crop = tuple(slice(int(v), int(v+size)) for v in lo)
            coverage = f['coverage'][crop] > 0
            support = f['support'][crop].astype('f8')
            rho = f['rho'][(slice(None), *crop)].astype('f8')
            assert np.isfinite(rho[:,coverage]).all()
            kernels = op.kernel((size,)*3, target_shape, target_lo-lo)
            shape = kernels.shape[1:]
            for mode in ('covered', 'fractional_support'):
                weight = coverage.astype('f8') if mode == 'covered' else support*coverage
                tensors = np.zeros((4, *target_shape, 6))
                spectra = [fft.rfftn(np.where(coverage, rho[d]-1, 0)*weight,
                                    s=shape, workers=16) for d in range(4)]
                for component, kernel in enumerate(kernels):
                    khat = fft.rfftn(kernel, workers=16)
                    for d in range(4):
                        tensors[d,...,component] = fft.irfftn(spectra[d]*khat, s=shape, workers=16)[:32,:32,:32]
                expected = (core_rho-1)*(1 if mode == 'covered' else support_core)
                trace_error = float(abs(tensors[...,[0,3,5]].sum(-1)-expected)[:,core_valid].max())
                assert trace_error < 1e-8*max(1., float(abs(expected[:,core_valid]).max()))
                eig = op.eigensystem(tensors)
                labels = (eig > .2).sum(-1).astype('u1')
                labels[:,~core_valid] = 255
                key = f'{mode}_{size}'
                layers[key] = (eig, labels)
                results.append(dict(key=key, source_width_mpc_h=size*CELL,
                    trace_maxabs=trace_error, supported_source_fraction=float((support>=.5).mean()),
                    fractions=[float((labels[:,core_valid]==c).mean()) for c in range(4)]))
                print(json.dumps(results[-1]), flush=True)
                del spectra
            del kernels, rho
        comparisons = []
        for mode in ('covered','fractional_support'):
            for left, right in ((96,144),(144,192)):
                a, la = layers[f'{mode}_{left}']; b, lb = layers[f'{mode}_{right}']
                comparisons.append(dict(left=f'{mode}_{left}', right=f'{mode}_{right}',
                    eigen_abs_quantiles=np.quantile(abs(a-b)[:,core_valid],[.5,.95,.99,1]).tolist(),
                    label_changed_fraction=float((la[:,core_valid]!=lb[:,core_valid]).mean())))
        a,la = layers['covered_192']; b,lb = layers['fractional_support_192']
        comparisons.append(dict(left='covered_192',right='fractional_support_192',
            eigen_abs_quantiles=np.quantile(abs(a-b)[:,core_valid],[.5,.95,.99,1]).tolist(),
            label_changed_fraction=float((la[:,core_valid]!=lb[:,core_valid]).mean())))
        np.savez_compressed(out/'tidal_sensitivity_layers.npz', rho=core_rho,
            origin_mpc_h=origin+target_lo*CELL, cell_mpc_h=CELL, support=support_core,valid=core_valid,
            **{k+'_eigen':v[0] for k,v in layers.items()},
            **{k+'_class':v[1] for k,v in layers.items()})
        rcrop=tuple(slice(int(v-44),int(v+44)) for v in center)
        rho=f['rho'][(slice(None),*rcrop)]
        valid=(f['support'][rcrop]>=.5)&(f['coverage'][rcrop]>0)
        pooled=rho.reshape(4,44,2,44,2,44,2).mean((2,4,6))
        pvalid=valid.reshape(44,2,44,2,44,2).all((1,3,5))
        matched_valid=np.repeat(np.repeat(np.repeat(pvalid,2,0),2,1),2,2)
        diagnostics=[]
        for d in range(4):
            peaks=(rho[d]==ndimage.maximum_filter(rho[d],size=3)) & valid & (rho[d]>2)
            interior=ndimage.binary_erosion(valid,iterations=2)
            diagnostics.append(dict(draw=d, density_quantiles=np.quantile(rho[d][valid],[.01,.5,.95,.99,1]).tolist(),
                native_std=float(rho[d][matched_valid].std()), pooled_std=float(pooled[d][pvalid].std()),
                local_maxima_above_2=int((peaks&interior).sum())))
        np.savez_compressed(out/'render_controls.npz',rho=rho,valid=valid,pooled=pooled,pvalid=pvalid)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, ListedColormap
        fig,axes=plt.subplots(2,3,figsize=(15,10),layout='constrained')
        images=[rho[0,:,:,44],np.repeat(np.repeat(pooled[0,:,:,22],2,0),2,1),rho.mean(0)[:,:,44]]
        titles=['Draw 1: native slice','Draw 1: 2-cell pooled slab','Mean of four: native slice']
        masks=[valid[:,:,44],matched_valid[:,:,44],valid[:,:,44]]
        for ax,im,title,mask in zip(axes[0],images,titles,masks):
            ax.imshow(np.ma.masked_where(~mask,im).T,origin='lower',norm=LogNorm(.5,8),cmap='magma',extent=[0,88*CELL,0,88*CELL],interpolation='nearest')
            ax.set_title(title); ax.set_xlabel('Mpc/h')
        for ax,key in zip(axes[1],('covered_96','covered_192','fractional_support_192')):
            ax.imshow(np.ma.masked_where(~core_valid[:,:,16],layers[key][1][0,:,:,16]).T,origin='lower',vmin=0,vmax=3,
                cmap=ListedColormap(['#17233b','#4c91ba','#eea34b','#ce4d83']),interpolation='nearest',extent=[0,32*CELL,0,32*CELL])
            ax.set_title(key+' / draw 1');ax.set_xlabel('Mpc/h')
        fig.suptitle('Resolution controls and provisional tidal sensitivity\nBottom: same central 216.5 Mpc/h region; void / wall / filament / knot, threshold 0.2; assumed mean-density exterior',fontsize=13)
        fig.savefig(out/'diagnostic.png',dpi=150);plt.close(fig)
    record=dict(schema='loa-structure-audit-v1',source=str(SOURCE),source_sha256=EXPECTED,
        script_sha256=sha(__file__),operator_sha256=sha(OPERATOR),job=os.environ['SLURM_JOB_ID'],
        numerical_checks=checks,contexts=results,comparisons=comparisons,morphology=diagnostics,
        center_cell=center.tolist(),target_shape=list(target_shape),valid_target_cells=int(core_valid.sum()),threshold=.2,
        support_audit=dict(total=core_valid.size,no_coverage=int((coverage_core==0).sum()),
            below_half=int((support_core<.5).sum()),zero_support=int((support_core==0).sum()),
            all_density_finite=bool(np.isfinite(core_rho).all())),
        status='numerically checked sensitivity layers; physical exterior and DESI calibration not validated',
        limitations=['No mock truth or calibrated real-DESI uncertainty',
            'Local maxima are density-grid peaks, not resolved halos or certified clusters',
            'Increasing finite source boxes is not a test of all unknown exterior mass',
            'Two support treatments are alternate assumptions, not posterior draws',
            'Pooled slice averages a two-cell slab; native slice is one cell',
            'No new smoothing, inference, training, or release replacement'])
    (out/'RESULTS.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    main(parser.parse_args().out)
