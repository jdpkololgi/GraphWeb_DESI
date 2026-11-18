"""Matched density/tidal cubes and Loa parent galaxies; no posterior sampling."""
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket

import h5py
import numpy as np
from scipy import fft

BASE=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1/production_v1')
ATLAS=Path('/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1')
OPERATOR=Path('/global/u2/d/dkololgi/TNG/Illustris/workflows/sbi/e2e_footprint_tidal_operator.py')
CELL=6.766
N=32

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

def main(out):
    assert os.environ.get('SLURM_JOB_ID') and socket.gethostname().startswith('nid')
    out.mkdir(parents=True,exist_ok=False)
    plan=json.loads((ATLAS/'PLAN.json').read_text())
    assert sha(ATLAS/'vac_properties.npz')==plan['properties_sha256']
    with np.load(ATLAS/'vac_properties.npz') as f:
        xyz,capids,ids,z=(f[k] for k in ('xyz','CAP','TARGETID','Z'))
    assert xyz.shape==(5436413,3) and np.isfinite(xyz).all()
    spec=importlib.util.spec_from_file_location('tidal',OPERATOR)
    op=importlib.util.module_from_spec(spec);spec.loader.exec_module(op)
    manifest=dict(schema='loa-cube-projections-v1',cell=CELL,shape=[N]*3,draws=4,
        model='Expanded U-Net · seed 17',contexts=[96,144],cubes=[],
        selection='Three anchor redshift ranges per cap. Choose a complete 32-cell cube near the coordinate-wise catalogue median; prefer support-rich geometry among nearest candidate anchors, without reading density values. Ranges select anchors, not galaxies within a cube.',
        galaxy_population='All Loa parent-VAC galaxies inside each half-open cube, 0.15 ≤ z < 0.55; no thinning, tile-edge or property cuts. These are observed context galaxies, not assigned T-web labels.',
        boundary='Free-space cubic-cell tidal operator, all prediction-covered source density, mean-density exterior and uncovered source cells. 96/144-cell source widths are sensitivity choices, not exterior posterior draws.',
        provenance=dict(parent_sha256=plan['properties_sha256'],operator_sha256=sha(OPERATOR),script_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID']))
    checks=[]
    def save(name,a):
        raw=np.ascontiguousarray(a).tobytes();packed=gzip.compress(raw,6,mtime=0)
        assert gzip.decompress(packed)==raw
        (out/name).write_bytes(packed)
        return dict(file=name,sha256=hashlib.sha256(packed).hexdigest(),bytes=len(packed))
    for cap,capid in [('NGC',1),('SGC',0)]:
        complete=json.loads((BASE/f'{cap}_seed17/COMPLETE.json').read_text())
        record=next(r for r in complete['results'] if r['model']=='unet_x_17')
        source=Path(record['field']['path']);assert complete['passed'] and sha(source)==record['field']['sha256']
        with h5py.File(source) as f:
            origin=np.asarray(f.attrs['origin_mpc_h'],dtype='f8');shape=np.array(f['coverage'].shape)
            assert f.attrs['cell_mpc_h']==CELL and f['rho'].shape[0]==4
            for j,(zlo,zhi) in enumerate(zip((.15,.28,.41),(.28,.41,.55))):
                rows=np.flatnonzero((capids==capid)&(z>=zlo)&(z<zhi))
                median=np.median(xyz[rows],axis=0)
                sample=np.random.default_rng(20261028+capid*10+j).choice(rows,min(4096,len(rows)),replace=False)
                candidates=sample[np.argsort(np.sum((xyz[sample]-median)**2,axis=1),kind='stable')[:128]]
                best=None;seen=set()
                for anchor in candidates:
                    center=np.floor((xyz[anchor]-origin)/CELL).astype(int)
                    lo=center-N//2;signature=tuple(lo)
                    if signature in seen or np.any(lo<0) or np.any(lo+N>shape):continue
                    seen.add(signature);crop=tuple(slice(int(a),int(a+N)) for a in lo)
                    if not np.all(f['coverage'][crop]>0):continue
                    support=f['support'][crop].astype('f8');score=float(support.mean())
                    if best is None or score>best[0]:best=(score,int(anchor),center,lo,crop,support)
                    if score>=.9:break
                assert best is not None,f'No fully computed cube for {cap} {zlo}'
                score,anchor,center,lo,crop,support=best
                rho=f['rho'][(slice(None),*crop)].astype('<f4')
                assert np.isfinite(rho).all() and np.all(rho>0)
                edge=origin+lo*CELL;upper=edge+N*CELL
                inside=(capids==capid)&np.all((xyz>=edge)&(xyz<upper),axis=1)
                selected=np.flatnonzero(inside)
                galaxies=np.column_stack((xyz[selected],z[selected])).astype('<f8')
                assert len(galaxies)>0
                cid=f'{cap}-{j+1}'
                item=dict(id=cid,label=f'{cap} · cube {j+1} · anchor z={z[anchor]:.3f}',cap=cap,
                    anchor_z=float(z[anchor]),anchor_redshift_range=[zlo,zhi],origin=edge.tolist(),
                    center=(edge+N*CELL/2).tolist(),shape=[N]*3,cell=CELL,galaxy_count=len(selected),
                    support_mean=score,source_sha256=record['field']['sha256'],files={})
                item['files']['rho']=save(cid+'-rho.bin',rho)
                item['files']['support']=save(cid+'-support.bin',support.astype('<f4'))
                item['files']['galaxies']=save(cid+'-galaxies.bin',galaxies)
                item['galaxy_targetids_sha256']=hashlib.sha256(ids[selected].tobytes()).hexdigest()
                eigenvalues={}
                for width in (96,144):
                    source_lo=np.maximum(0,center-width//2);source_hi=np.minimum(shape,center+width//2)
                    scrop=tuple(slice(int(a),int(b)) for a,b in zip(source_lo,source_hi))
                    covered=f['coverage'][scrop]>0
                    source_rho=f['rho'][(slice(None),*scrop)].astype('f8')
                    assert np.isfinite(source_rho[:,covered]).all()
                    kernels=op.kernel(tuple(source_hi-source_lo),(N,)*3,lo-source_lo)
                    pad=kernels.shape[1:]
                    tensor=np.empty((4,N,N,N,6),dtype='f8')
                    spectra=[fft.rfftn(np.where(covered,source_rho[d]-1,0),s=pad,workers=16) for d in range(4)]
                    for component,k in enumerate(kernels):
                        khat=fft.rfftn(k,workers=16)
                        for d in range(4):tensor[d,...,component]=fft.irfftn(spectra[d]*khat,s=pad,workers=16)[:N,:N,:N]
                    trace=float(np.max(abs(tensor[...,[0,3,5]].sum(-1)-(rho.astype('f8')-1))))
                    assert trace<1e-8
                    eig=op.eigensystem(tensor);e32=eig.astype('<f4')
                    assert np.all(np.diff(e32,axis=-1)>=0) and np.isfinite(e32).all()
                    error=float(np.max(abs(eig-e32)));assert error<1e-6
                    eigenvalues[width]=e32
                    item['files'][f'eigen_{width}']=save(f'{cid}-eigen-{width}.bin',e32)
                    checks.append(dict(cube=cid,width=width,trace_maxabs=trace,float32_maxabs=error,
                        source_shape=(source_hi-source_lo).tolist(),source_covered_fraction=float(covered.mean())))
                    del spectra,kernels,tensor,source_rho
                for eig in eigenvalues.values():
                    for threshold in (0.,.2,.5,1.):
                        labels=(eig.astype('f8')>threshold).sum(-1)
                        for axis in range(1,4):
                            fractions=np.stack([(labels==c).mean(axis=axis) for c in range(4)])
                            assert np.max(abs(fractions.sum(0)-1))<1e-12
                item['context_label_change_at_02']=float(((eigenvalues[96]>.2).sum(-1)!=(eigenvalues[144]>.2).sum(-1)).mean())
                # Independent reference pixels for browser projection checks.
                labels=(eigenvalues[144][0].astype('f8')>.2).sum(-1)
                item['reference']=dict(axis=2,start=12,stop=20,draw=0,threshold=.2,
                    mean_rho=float(rho[0,10,11,12:20].astype('f8').mean()),
                    class_fractions=[float((labels[10,11,12:20]==c).mean()) for c in range(4)])
                manifest['cubes'].append(item)
                print(json.dumps(dict(cube=cid,galaxies=len(selected),support_mean=score,context_changed=item['context_label_change_at_02'])),flush=True)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False))
    (out/'QA.json').write_text(json.dumps(dict(passed=True,checks=checks,job=os.environ['SLURM_JOB_ID']),indent=2))
    print(json.dumps(dict(passed=True,cubes=len(manifest['cubes']),bytes=sum(p.stat().st_size for p in out.iterdir()))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args();main(args.out)
