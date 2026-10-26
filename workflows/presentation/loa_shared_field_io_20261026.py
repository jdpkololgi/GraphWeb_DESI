"""Bounded shared-density storage and independent checks; no model fitting."""
import hashlib
import itertools
import json
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np

CELL = 6.766


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path, value, replace=False):
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    with os.fdopen(fd, 'w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.flush()
        os.fsync(stream.fileno())
    if replace:
        os.replace(temporary, path)
    else:
        os.link(temporary, path)
        os.unlink(temporary)


def coverage_for(geometry):
    if len(geometry.starts) >= np.iinfo(np.uint16).max:
        raise ValueError('Coverage counter overflow')
    coverage = np.zeros(geometry.shape, dtype='u2')
    for i in range(len(geometry.starts)):
        coverage[geometry.crop(i)[2:]] += 1
    return coverage


def blocks(values):
    batch, nx, ny, nz = values.shape
    return values.reshape(batch, nx // 4, 4, ny // 4, 4, nz // 4, 4)


def check_coverage(coverage, points, origin, expected_rows):
    """Every parent galaxy, not only independently cropped tile interiors."""
    shape = np.asarray(coverage.shape)
    if np.any(shape % 4) or len(points) != expected_rows:
        raise ValueError('Coverage shape or parent-row count mismatch')
    valid = coverage > 0
    reduced = blocks(valid[None])
    if not np.array_equal(reduced.min((2, 4, 6)), reduced.max((2, 4, 6))):
        raise ValueError('Partial mass blocks in coverage')
    for start in range(0, len(points), 100000):
        xyz = np.asarray(points[start:start + 100000])
        if not np.isfinite(xyz).all():
            raise ValueError('Nonfinite catalogue coordinates')
        lower = np.floor((xyz - origin) / CELL - .5).astype(np.int64)
        if np.any(lower < 0) or np.any(lower + 1 >= shape):
            raise ValueError('Parent galaxy outside interpolation bounds')
        for offset in itertools.product((0, 1), repeat=3):
            cell = lower + offset
            if not np.all(valid[cell[:, 0], cell[:, 1], cell[:, 2]]):
                raise ValueError('Parent galaxy has an uncovered interpolation corner')
    return dict(parent_rows=len(points), all_eight_corners_covered=True,
                covered_native_cells=int(valid.sum()))


def merge_support(destination, seen, value, crop):
    """Support is the unchanged observed random-support fraction, not density."""
    value = np.asarray(value, dtype='f4')
    if value.shape != destination[crop].shape or not np.isfinite(value).all():
        raise ValueError('Invalid support shape or values')
    if np.any(value < 0) or np.any(value > 1):
        raise ValueError('Support fraction outside [0,1]')
    shared = seen[crop]
    if np.any(shared) and not np.array_equal(destination[crop][shared], value[shared]):
        raise ValueError('Inconsistent support in overlapping observation windows')
    destination[crop] = value
    seen[crop] = True


def write_field(path, z, u, coverage, support, coarse_geometry, fine_geometry,
                chart, grid_origin, cap, model):
    """Decode block-aligned slabs, then verify float32 storage before publication."""
    import torch
    from workflows.sbi.e2e_vdm_context_models import decode_density, block_mean

    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    temporary = path.with_name('.' + path.name + '.partial')
    fg, cg = fine_geometry, coarse_geometry
    offset = np.asarray(fg.origin) // 4 - np.asarray(cg.origin)
    if np.any(np.asarray(fg.origin) % 4) or np.any(offset < 0):
        raise ValueError('Invalid coarse/fine physical alignment')
    if np.any(offset + np.asarray(fg.shape) // 4 > cg.shape):
        raise ValueError('Fine canvas exceeds coarse canvas')
    if tuple(u.shape) != (len(z), 1, *fg.shape) or tuple(z.shape[2:]) != cg.shape:
        raise ValueError('Latent shape mismatch')
    if coverage.shape != fg.shape or support.shape != fg.shape:
        raise ValueError('Mask shape mismatch')
    if not np.isfinite(support).all() or np.any((support < 0) | (support > 1)):
        raise ValueError('Invalid support')
    if np.any(support[coverage == 0] != 0):
        raise ValueError('Support outside covered windows')
    origin = np.asarray(grid_origin) + np.asarray(fg.origin) * CELL
    mass = 0.
    with h5py.File(temporary, 'x') as handle:
        rho_out = handle.create_dataset('rho', shape=(len(z), *fg.shape), dtype='f4',
            chunks=(1, *(min(n, k) for n, k in zip(fg.shape, (16, 32, 32)))), compression='lzf')
        handle.create_dataset('coverage', data=coverage, compression='lzf')
        handle.create_dataset('support', data=support, compression='lzf')
        handle.create_dataset('coarse_logrho_standardized', data=z[:, 0].cpu().numpy(), compression='lzf')
        handle.attrs.update(origin_mpc_h=origin, cell_mpc_h=CELL, distance_unit='Mpc/h',
            cap=cap, model=model, draws=len(z), nfe=128, shared_state=True,
            posterior_qualified=False, coarse_origin_raw_cells=np.asarray(cg.origin) * 8,
            fine_origin_raw_cells=np.asarray(fg.origin) * 2,
            support_semantics='unchanged pooled random-support fraction; zero is not filled')
        for lo in range(0, fg.shape[0], 16):
            hi = min(lo + 16, fg.shape[0])
            begin = offset + [lo // 4, 0, 0]
            size = np.array([hi - lo, fg.shape[1], fg.shape[2]]) // 4
            crop = (slice(None), slice(None), *(slice(int(a), int(a + n)) for a, n in zip(begin, size)))
            coarse = torch.exp(z[crop].double() * chart['coarse_logrho']['std'][0]
                               + chart['coarse_logrho']['mean'][0])
            rho = decode_density(coarse, u[:, :, lo:hi].double() * chart['fine_residual']['std'][0])
            if not torch.isfinite(rho).all() or not (rho > 0).all():
                raise ValueError('Nonpositive or nonfinite decoded density')
            mass = max(mass, float(((block_mean(rho) - coarse) / coarse).abs().max()))
            values = rho[:, 0].cpu().numpy().astype('f4')
            values[:, coverage[lo:hi] == 0] = np.nan
            rho_out[:, lo:hi] = values
        handle.flush()
    with temporary.open('rb') as stream:
        os.fsync(stream.fileno())
    if mass >= 2e-6:
        raise ValueError('Decoded mass gate failed')
    stored = verify_field(temporary, chart, origin, offset, coverage, len(z))
    # Hard-link publication is atomic and refuses to replace any existing field.
    os.link(temporary, path)
    os.unlink(temporary)
    return dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size,
                decoded_mass_max_relative_error=mass, **stored)


def verify_field(path, chart, expected_origin, expected_offset, expected_coverage, draws):
    """Independent NumPy verification of persisted HDF5, including rounding."""
    maximum = 0.
    with h5py.File(path, 'r') as handle:
        coverage = handle['coverage'][:]
        shape = np.asarray(coverage.shape)
        if not np.array_equal(coverage, expected_coverage):
            raise ValueError('Stored coverage mismatch')
        if (not np.array_equal(handle.attrs['origin_mpc_h'], expected_origin)
                or handle.attrs['cell_mpc_h'] != CELL or not handle.attrs['shared_state']
                or handle.attrs['posterior_qualified'] or handle['rho'].shape != (draws, *shape)):
            raise ValueError('Stored field metadata mismatch')
        offset = np.asarray(expected_offset, dtype=int)
        for lo in range(0, shape[0], 16):
            hi = min(lo + 16, shape[0])
            rho = handle['rho'][:, lo:hi].astype('f8')
            valid = coverage[lo:hi] > 0
            if not np.array_equal(np.isfinite(rho), np.broadcast_to(valid, rho.shape)):
                raise ValueError('Stored finite/coverage mask mismatch')
            if not np.all(rho[:, valid] > 0):
                raise ValueError('Stored density is not strictly positive')
            support = handle['support'][lo:hi]
            if (not np.isfinite(support).all() or np.any((support < 0) | (support > 1))
                    or np.any(support[~valid] != 0)):
                raise ValueError('Stored survey support mismatch')
            block_valid = blocks(valid[None]).min((2, 4, 6))[0]
            if not np.array_equal(block_valid, blocks(valid[None]).max((2, 4, 6))[0]):
                raise ValueError('Stored mask splits a coarse mass block')
            z = handle['coarse_logrho_standardized'][:,
                offset[0] + lo // 4:offset[0] + hi // 4,
                offset[1]:offset[1] + shape[1] // 4,
                offset[2]:offset[2] + shape[2] // 4].astype('f8')
            coarse = np.exp(z * chart['coarse_logrho']['std'][0] + chart['coarse_logrho']['mean'][0])
            means = blocks(rho).mean((2, 4, 6))
            if np.any(block_valid):
                error = np.abs(means[:, block_valid] - coarse[:, block_valid]) / coarse[:, block_valid]
                if not np.isfinite(error).all():
                    raise ValueError('Nonfinite stored mass error')
                maximum = max(maximum, float(error.max()))
    if maximum >= 2e-6:
        raise ValueError('Stored float32 mass gate failed')
    return dict(stored_float32_mass_max_relative_error=maximum, stored_field_verified=True)
