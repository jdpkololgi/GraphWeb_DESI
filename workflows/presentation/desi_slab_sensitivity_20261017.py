#!/usr/bin/env python3
"""Preliminary nested-slab check from saved plotted rows, not a 3D closure test."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json, hashlib
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010/filament_and_local_geometry_20261017')
p=ROOT/'data/filament_audit_rows.npz'
d=dict(np.load(p)); xy=np.c_[d['x_Mpc'],d['y_Mpc']]; h=d['height_Mpc']; radius=7/.6766
interior=d['projected_shape_valid'] & (abs(h)<2.5)
cut=float(json.loads((ROOT/'data/filament_audit.json').read_text())['local_shape']['top20pct_cut'])
rows=[]; arrays={}
for depth in (5,10,20):
 mask=abs(h)<depth/2; ix=np.flatnonzero(mask); tree=cKDTree(xy[ix]); shape=np.full(len(h),np.nan); count=np.zeros(len(h),int)
 for i,near in zip(ix,tree.query_ball_point(xy[ix],radius)):
  q=ix[near]; q=q[q!=i]; count[i]=len(q)
  if len(q)>=8:
   e=np.linalg.eigvalsh(np.cov(xy[q].T));shape[i]=(e[1]-e[0])/sum(e)
 valid=interior & np.isfinite(shape); linear=valid & (shape>=cut); low=d['P_WEB'][:,2]<.2
 rows.append(dict(depth_mpc=depth,displayed_n=int(mask.sum()),fixed_target_pool_n=int(interior.sum()),valid_fixed_targets=int(valid.sum()),line_like_fixed_cut=int(linear.sum()),low_pfil_line_like=int((linear&low).sum()),low_fraction_line_like=float(low[linear].mean()),median_neighbors_fixed_targets=float(np.median(count[interior]))))
 arrays[f'shape_depth_{depth}']=shape;arrays[f'neighbors_depth_{depth}']=count
# Paired comparison on targets with sufficient neighbours at every depth.
common=interior.copy()
for depth in (5,10,20):common &= np.isfinite(arrays[f'shape_depth_{depth}'])
paired=[]
for depth in (5,10,20):
 line=common & (arrays[f'shape_depth_{depth}']>=cut)
 paired.append(dict(depth_mpc=depth,common_valid_targets=int(common.sum()),line_like=int(line.sum()),low_pfil=int((line & (d['P_WEB'][:,2]<.2)).sum()),low_fraction=float((d['P_WEB'][line,2]<.2).mean())))
result=dict(source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),fixed_shape_threshold=cut,fixed_targets='Original valid shape targets with |height|<2.5 Mpc; x/y edge margin retained',radius_mpc=radius,nested_slab_results=rows,paired_common_valid=paired,full_3d_spheres_available=0,why_no_full_spheres='Original half-depth 10 Mpc is smaller than diagnostic radius 10.34585 Mpc; no original-slab target has a complete sphere.',conclusion='Checks sensitivity within the existing slab only. Cannot establish or exclude a slab-selection explanation without off-slab neighbours and held-out mock truth. No retraining decision justified by this statistic alone.')
(ROOT/'data/slab_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n')
np.savez_compressed(ROOT/'data/slab_sensitivity_arrays.npz',TARGETID=d['TARGETID'],common_valid=common,**arrays)
print(json.dumps(result,indent=2))
