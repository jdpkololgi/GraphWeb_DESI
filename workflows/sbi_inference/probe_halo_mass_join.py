#!/usr/bin/env python3
"""Feasibility probe for joining halo mass onto path1 mock galaxies via
(FILE_NUM, BOX_INDEX) -> CompaSO halo_info, to anchor 'cluster' to mass rather than
the smoothing-dependent lambda1 threshold."""
from __future__ import annotations
import numpy as np
import fitsio

WEDGE = ("/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/"
         "path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3_wedge_targets.fits")
HALO_DIR = "/pscratch/sd/d/dkololgi/AbacusSummit_densities/AbacusSummit_base_c000_ph000/halos/z0.200/halo_info"

t = fitsio.read(WEDGE, columns=["FILE_NUM", "BOX_INDEX", "HALO_INDEX"])
fn, bi, hi = t["FILE_NUM"], t["BOX_INDEX"], t["HALO_INDEX"]
print(f"path1 wedge N={len(t)}")
for nm, a in [("FILE_NUM", fn), ("BOX_INDEX", bi), ("HALO_INDEX", hi)]:
    print(f"  {nm}: min={a.min()} max={a.max()} frac<0={np.mean(a<0):.3f} nunique={len(np.unique(a))}")

# load one CompaSO slab, get N (particle count) + particle mass from header
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
cat = CompaSOHaloCatalog(f"{HALO_DIR}/halo_info_000.asdf", fields=["N"], cleaned=False)
N = np.asarray(cat.halos["N"])
pm = cat.header.get("ParticleMassHMsun", cat.header.get("ParticleMassMsun", None))
print(f"\nhalo_info_000: nhalos={len(N)}  N range [{N.min()},{N.max()}]  ParticleMass={pm:.3e}")
print(f"  -> mass range [{N.min()*pm:.2e}, {N.max()*pm:.2e}] Msun/h")
print(f"  nhalos with M>1e13 (cluster-ish): {(N*pm>1e13).sum()}")

# does (FILE_NUM=0, BOX_INDEX) index into slab 0?
m0 = fn == 0
print(f"\npath1 FILE_NUM==0: {m0.sum()} gals; BOX_INDEX max={bi[m0].max()} vs slab0 nhalos={len(N)} "
      f"-> in range: {bi[m0].max() < len(N)}")
print(f"  HALO_INDEX (FILE_NUM==0) max={hi[m0].max()}  (which one indexes the slab?)")
# sample masses via BOX_INDEX
valid = m0 & (bi >= 0) & (bi < len(N))
if valid.sum() > 100:
    mass = N[bi[valid]] * pm
    print(f"  via BOX_INDEX: galaxy host-halo mass median={np.median(mass):.2e}, "
          f"frac>1e13={np.mean(mass>1e13):.3f}, max={mass.max():.2e}")
