#!/usr/bin/env python3
"""Validate the (FILE_NUM, BOX_INDEX/HALO_INDEX) -> CompaSO halo indexing by checking
that the joined halo position matches the authoritative host_halos_unique x_com
(built by the user's pipeline with key=(FILE_NUM, BOX_INDEX)). Tests index column x
cleaned flag to find the convention that reproduces x_com."""
from __future__ import annotations
import numpy as np
import fitsio

PRE = ("/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/"
       "path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3")
HALO_DIR = "/pscratch/sd/d/dkololgi/AbacusSummit_densities/AbacusSummit_base_c000_ph000/halos/z0.200/halo_info"
HK = "/pscratch/sd/d/dkololgi/abacus/host_halos_unique_keys.npy"
HX = "/pscratch/sd/d/dkololgi/abacus/host_halos_unique_xcom_points.npy"

tg = fitsio.read(f"{PRE}_wedge_targets.fits", columns=["FILE_NUM", "BOX_INDEX", "HALO_INDEX"])
keys = np.load(HK)            # [M,2] (FILE_NUM, BOX_INDEX)
xcom = np.load(HX, mmap_mode="r")  # [M,3]
print("host_halos keys", keys.shape, "xcom", xcom.shape)

f0 = 21
gsel = tg["FILE_NUM"] == f0
print(f"path1 FILE_NUM=={f0}: {gsel.sum()} galaxies")

# authoritative key->xcom for slab f0
ksel = keys[:, 0] == f0
kb = keys[ksel, 1]; kx = np.asarray(xcom[ksel])
kmap = {int(b): i for i, b in enumerate(kb)}
print(f"host_halos slab {f0}: {ksel.sum()} unique host halos; BOX_INDEX range [{kb.min()},{kb.max()}]")

from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
samp = np.where(gsel)[0][:2000]
for cleaned in (False, True):
    try:
        cat = CompaSOHaloCatalog(f"{HALO_DIR}/halo_info_{f0:03d}.asdf",
                                 fields=["N", "x_L2com"], cleaned=cleaned)
    except Exception as e:
        print(f"cleaned={cleaned}: load failed: {e}"); continue
    cols = list(cat.halos.colnames)
    xL = np.asarray(cat.halos["x_L2com"]); Ncat = np.asarray(cat.halos["N"])
    print(f"\ncleaned={cleaned}: nhalos={len(xL)} cols~{[c for c in cols if 'x_' in c or c=='N'][:6]}")
    for idxcol in ("BOX_INDEX", "HALO_INDEX"):
        bvals = tg[idxcol][samp]
        ok = bvals < len(xL)
        # expected xcom from host_halos for these galaxies (via their BOX_INDEX key)
        exp = np.array([kmap.get(int(b), -1) for b in tg["BOX_INDEX"][samp]])
        m = ok & (exp >= 0)
        if m.sum() < 50:
            print(f"  {idxcol}: too few matchable ({m.sum()})"); continue
        d = np.linalg.norm(xL[bvals[m]] - kx[exp[m]], axis=1)
        print(f"  index={idxcol}: median |x_compaso - x_host| = {np.median(d):.3f} Mpc/h "
              f"(match if ~0)  [frac<1Mpc/h: {np.mean(d<1):.2f}]")
