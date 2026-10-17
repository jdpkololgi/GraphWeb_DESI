"""Independent result reconciliation and secondary paired summaries."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from prepare_parent import digest

def main():
 p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);a=p.parse_args();r=json.loads((a.results/'RESULTS.json').read_text());g=dict(np.load(a.results/'geometry.npz'));assert digest(a.results/'geometry.npz')==r['provenance']['geometry_sha256']
 checks=0
 for row in r['rows']:
  k=int(np.argmin(abs(np.array([.5,1,1.5])*7/.6766-row['radius_Mpc'])));c=g['cohort_'+row['cohort']];valid=c&g[f'eligible_{k}']&(g[f'nfull_{k}']>=row['min_neighbors']);v=g[f'full_{k}'][valid]
  assert len(v)==row['full']['n'];assert np.all(np.isfinite(v));np.testing.assert_allclose(v.sum(1),1,atol=1e-12);np.testing.assert_array_equal(np.bincount(v.argmax(1),minlength=3),row['full'].get('dominant_counts',[0,0,0]));checks+=4
 for k in range(3):assert np.all(g[f'nslab_{k}']<=g[f'nfull_{k}']);checks+=1
 extra={}
 for name in ['low_fil_wall','low_fil_void','line_fil_argmax']:
  c=g['cohort_'+name];v=c&g['eligible_1']&(g['nslab_1']>=8);s=g['slab_1'][v].argmax(1);f=g['full_1'][v].argmax(1);transition=np.zeros((3,3),int);np.add.at(transition,(s,f),1)
  common=c.copy()
  for depth in [20,40,80]:common&=np.isfinite(g[f'shape_a0_o0_d{depth}'])
  extra[name]={'paired_3d_shape_transition_slab_rows_full_columns_LPS':transition.tolist(),'same_targets_depth20_40_80_n':int(common.sum()),'same_targets_depth20_40_80_line_counts':[int(np.sum(common&(g[f'shape_a0_o0_d{depth}']>=.5847730874520835))) for depth in [20,40,80]],'median_full_neighbors_main':float(np.median(g['nfull_1'][c]))}
 # Normal direction is descriptive; low |normal_z| means the sheet intersects the original slab edge-on.
 c=g['cohort_low_fil_wall']&g['eligible_1']&(g['nfull_1']>=8)&(g['full_1'].argmax(1)==1)
 extra['sheet_normal']={'n':int(c.sum()),'median_absolute_dot_original_slab_normal':float(np.median(abs(g['normal_1'][c,2]))),'fraction_abs_dot_lt_half':float(np.mean(abs(g['normal_1'][c,2])<.5))}
 (a.results/'VALIDATION.json').write_text(json.dumps({'checks':checks,'secondary_summaries':extra,'note':'20/40/80 paired diagnostic is secondary after all-five-depth common cohort proved tiny. No threshold or target definition changed.'},indent=2,allow_nan=False)+'\n');print(json.dumps(extra))
if __name__=='__main__':main()
