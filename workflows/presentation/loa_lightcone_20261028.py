"""Export continuous native Loa density in spherical slices, without new inference.

Run only on an approved compute node. One-degree declination quadrature and
half-degree RA sampling are display resampling, not improved field resolution.
"""
import argparse
import gzip
import hashlib
import json
import os
import socket
from pathlib import Path

import fitsio
import h5py
import numpy as np

BASE = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1/production_v1')
PARENT = Path('/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1/vac_properties.npz')
CAT = Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def north(ra, dec):
    # IAU north Galactic pole, consistent with the existing cap convention.
    a, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.sin(d)*np.sin(np.deg2rad(27.12825)) + np.cos(d)*np.cos(np.deg2rad(27.12825))*np.cos(a-np.deg2rad(192.85948)) > 0


def main(out):
    assert os.environ.get('SLURM_JOB_ID') and socket.gethostname().startswith('nid')
    out.mkdir(parents=True, exist_ok=False)
    assets = {}
    def save(name, array):
        raw = np.ascontiguousarray(array).tobytes()
        packed = gzip.compress(raw, 5, mtime=0)
        assert gzip.decompress(packed) == raw
        (out/name).write_bytes(packed)
        assets[name] = dict(bytes=len(packed), sha256=hashlib.sha256(packed).hexdigest())
        return name

    plan = json.loads(PARENT.with_name('PLAN.json').read_text())
    assert sha(PARENT) == plan['properties_sha256']
    with np.load(PARENT) as p:
        z = p['Z']; xyz = p['xyz']; radius = np.linalg.norm(xyz, axis=1)
        ra, dec = p['RA'], p['DEC']
        unit = np.column_stack((np.cos(np.deg2rad(dec))*np.cos(np.deg2rad(ra)),
                                np.cos(np.deg2rad(dec))*np.sin(np.deg2rad(ra)), np.sin(np.deg2rad(dec))))
        assert np.max(abs(xyz-radius[:, None]*unit)) < 1e-6
        assert np.all(north(ra, dec) == (p['CAP'] == 1))
        # Preserve the actual inference distance convention rather than guessing cosmology.
        order = np.argsort(z); zs, rs = z[order], radius[order]
        unique = np.r_[True, np.diff(zs) > 0]; zs, rs = zs[unique], rs[unique]
        assert np.all(np.diff(rs) >= -1e-7)
        knots = np.unique(np.r_[np.arange(0, len(zs), 100), len(zs)-1])
        zz, rr = zs[knots], rs[knots]
        err = float(np.max(abs(np.interp(zs, zz, rr)-rs)))
        assert err < .01, err
        # Extend only the sub-sample endpoint gaps to the declared selection boundaries.
        lo = rr[0] + (.15-zz[0])*(rr[1]-rr[0])/(zz[1]-zz[0])
        hi = rr[-1] + (.55-zz[-1])*(rr[-1]-rr[-2])/(zz[-1]-zz[-2])
        zz, rr = np.r_[.15, zz, .55], np.r_[lo, rr, hi]
        del xyz, unit, z, ra, dec, radius, order, zs, rs
    nr = int(np.ceil((hi-lo)/6.766))
    re = np.linspace(lo, hi, nr+1); rc = (re[1:]+re[:-1])/2
    ae = np.linspace(0, 360, 721); ac = (ae[1:]+ae[:-1])/2
    de = np.arange(-90, 91); dc = (de[1:]+de[:-1])/2
    # Complete observed bright catalogue, not the thinned homepage overlay.
    chunks = []; input_rows = 0
    stat = CAT.stat()
    with fitsio.FITS(CAT) as f:
        for start in range(0, f[1].get_nrows(), 250000):
            d = f[1].read(rows=np.arange(start, min(start+250000, f[1].get_nrows())),
                columns=['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE','BGS_TARGET'])
            input_rows += len(d)
            good = ((d['BGS_TARGET'] & 2) != 0) & (d['ZWARN'] == 0) & (d['DELTACHI2'] >= 25)
            good &= np.char.strip(d['SPECTYPE'].astype('U')) == 'GALAXY'
            good &= np.isfinite(d['RA']) & np.isfinite(d['DEC']) & (d['Z_not4clus'] >= .15) & (d['Z_not4clus'] < .55)
            chunks.append(d[good])
    gal = np.concatenate(chunks); ids = gal['TARGETID']
    assert len(np.unique(ids)) == len(ids), 'Duplicate TARGETIDs must be resolved explicitly'
    assert (CAT.stat().st_size, CAT.stat().st_mtime_ns) == (stat.st_size, stat.st_mtime_ns)
    cap = north(gal['RA'], gal['DEC'])
    manifest = dict(schema='loa-lightcone-v1', model='Expanded U-Net · seed 17', draws=4,
        z_range=[.15,.55], radial_edges=re.tolist(), radial_redshifts=np.interp(re, rr, zz).tolist(),
        ra_edges=ae.tolist(), dec_edges=de.tolist(), native_cell=6.766, caps={}, assets=assets,
        distance=dict(method='Monotonic interpolation of the actual parent inference coordinates', max_error_mpc_h=err),
        selection='Loa v2.1 BGS_BRIGHT_full_HPmapcut: BGS_TARGET & 2, ZWARN=0, DELTACHI2>=25, SPECTYPE=GALAXY, 0.15<=Z_not4clus<0.55. Unique TARGETIDs; no display thinning, magnitude recut, or field-support cut on galaxies.',
        projection='At each RA/radius sample, native cell lookup at one-degree declination-bin midpoints. Combine selected bins with exact solid-angle weights sin(dec_hi)-sin(dec_lo). Blank if any selected bin lacks positive random support or prediction coverage. No interpolation across missing regions, no added smoothing. This is an angular-thickness average, not a Cartesian slab or a connectivity measurement.',
        support='Positive native pooled random-support fraction AND prediction coverage. Random support is a finite-resolution footprint approximation, not a calibrated reliability probability. Coverage fraction is available separately; partially covered rays are blank by default.',
        provenance=dict(parent_sha256=plan['properties_sha256'], catalogue=dict(path=str(CAT),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns),
                        script_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID']))
    checks=[]
    for name, cap_value in [('NGC',True),('SGC',False)]:
        g = gal[cap == cap_value]; gi = np.floor(g['DEC']+90).astype(int)
        assert np.all((gi >= 0) & (gi < 180))
        rec = json.loads((BASE/f'{name}_seed17/COMPLETE.json').read_text())
        entry = next(x for x in rec['results'] if x['model']=='unet_x_17')
        source = Path(entry['field']['path']); assert rec['passed'] and sha(source)==entry['field']['sha256']
        slices=[]
        with h5py.File(source) as f:
            assert f.attrs['shared_state'] and f.attrs['cell_mpc_h']==6.766
            origin=np.asarray(f.attrs['origin_mpc_h']); shape=np.array(f['coverage'].shape)
            coverage=f['coverage'][:]; support=f['support'][:]; rho=f['rho'][:]
            assert rho.shape==(4,*shape)
            for k,d in enumerate(dc):
                a=np.deg2rad(ac)[None,:]; r=rc[:,None]; dd=np.deg2rad(d)
                x=np.broadcast_arrays(r*np.cos(dd)*np.cos(a),r*np.cos(dd)*np.sin(a),r*np.sin(dd))
                ix=np.floor((np.stack(x,axis=-1)-origin)/6.766).astype(int)
                valid=np.all((ix>=0)&(ix<shape),axis=-1) & (north(ac,d)[None,:]==cap_value)
                clipped=np.clip(ix,0,shape-1); loc=tuple(clipped[...,j] for j in range(3))
                valid &= (coverage[loc]>0) & (support[loc]>0)
                values=np.where(valid[None,...],rho[(slice(None),*loc)],0).astype('<f4')
                assert np.isfinite(values).all() and np.all(values[:,valid]>0)
                rows=g[gi==k]
                if not valid.any() and not len(rows): continue
                ga=np.column_stack((rows['RA'],rows['DEC'],np.interp(rows['Z_not4clus'],zz,rr),rows['Z_not4clus'])).astype('<f8')
                item=dict(bin=k,dec=[int(de[k]),int(de[k+1])],galaxies=len(rows),valid_pixels=int(valid.sum()),
                    density=save(f'{name}-{k}-rho.bin',values),galaxy_file=save(f'{name}-{k}-galaxies.bin',ga),
                    targetids_sha256=hashlib.sha256(rows['TARGETID'].tobytes()).hexdigest())
                if valid.any():
                    v,u=np.argwhere(valid)[len(np.argwhere(valid))//2]
                    item['reference']=dict(radial=int(v),ra=int(u),rho=values[:,v,u].tolist(),native_index=ix[v,u].tolist())
                slices.append(item)
            del coverage,support,rho
        assert sum(s['galaxies'] for s in slices)==len(g)
        manifest['caps'][name]=dict(slices=slices,galaxies=len(g),field_sha256=entry['field']['sha256'])
        checks.append(dict(cap=name,galaxies=len(g),slices=len(slices),reference_pixels=True))
        print(json.dumps(checks[-1]),flush=True)
    assert sum(x['galaxies'] for x in manifest['caps'].values())==len(gal)
    qa=dict(passed=True,checks=checks,input_catalogue_rows=input_rows,unique_bright_galaxies=len(gal),distance_error_mpc_h=err,
            job=os.environ['SLURM_JOB_ID'],byte_roundtrips=True,source_hashes=True)
    (out/'manifest.json').write_text(json.dumps(manifest,separators=(',',':')))
    (out/'QA.json').write_text(json.dumps(qa,indent=2))
    print(json.dumps(qa),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);main(p.parse_args().output)
