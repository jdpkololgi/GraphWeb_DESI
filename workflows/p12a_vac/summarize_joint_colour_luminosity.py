"""Use identical z/r cells for both parent variants when comparing colour residuals."""
import json,hashlib
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[2];root=R/'docs/evidence/p12a_joint_colour_luminosity_20260928'
h=np.load(root/'HISTOGRAMS.npz');x=json.loads((root/'RESULTS.json').read_text());rows=[];cen=(h['colour_edges'][:-1]+h['colour_edges'][1:])/2
for i,p in enumerate(['S','N']):
 for j in range(4):
  hs={n:h[n][i,j*10:(j+1)*10] for n in ['loa','v0.1','v1']};counts={n:a.sum(-1) for n,a in hs.items()};ok=np.logical_and.reduce([a>=20 for a in counts.values()]);a=hs['loa'][ok].sum(0);a/=a.sum()
  for v in ['v0.1','v1']:
   w=np.divide(counts['loa'],counts[v],out=np.zeros_like(counts['loa']),where=ok);b=(hs[v]*w[...,None]).sum((0,1));b/=b.sum()
   rows.append(dict(mock=v,photsys=p,zlo=float(h['z_edges'][j*10]),zhi=float(h['z_edges'][(j+1)*10]),matched_cells=int(ok.sum()),loa_weight_fraction_retained=float(counts['loa'][ok].sum()/counts['loa'].sum()),mean_colour_mock_minus_loa=float((b-a)@cen),colour_total_variation=float(.5*np.abs(b-a).sum())))
x['shared_rows']=rows;x['shared_support']='Identical cells with >=20 weighted/sample counts in all three catalogues; paired rows also retained.';x['summary_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(root/'RESULTS.json').write_text(json.dumps(x,indent=2)+'\n')
for r in rows:
 if r['zlo']>.4:print(r)
