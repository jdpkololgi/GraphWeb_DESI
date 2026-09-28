"""Compare pinned cumulative generator LFs with published differential r LFs.
Small table-only diagnostic; never alters catalogues or fits model parameters.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from astropy.table import Table
R=Path(__file__).resolve().parents[2]
O=R/'docs/evidence/p12a_published_lf_20260928'
rows=[]; sensitivity=[]; hashes={}
for branch in ['wsys','nowsys']:
 p=R/f'docs/evidence/p12a_photometry_swap_20260928/lf_{branch}.dat'
 x,y=np.loadtxt(p,unpack=True); assert np.all(np.diff(x)>0)
 assert np.all(np.diff(y)>=0)
 def cumulative(m):
  assert np.all((np.asarray(m)>=x[0])&(np.asarray(m)<=x[-1]))
  return 10**np.interp(m,x,y)
 hashes[str(p.relative_to(R))]=hashlib.sha256(p.read_bytes()).hexdigest()
 for f in sorted(O.glob('*.fits')):
  t=Table.read(f);m=np.array(t['mag']);width=np.diff(m)
  assert np.allclose(width,.25)
  model=(cumulative(m+.125)-cumulative(m-.125))/.25
  for lo,hi in [(-21,-20),(-22,-21),(-23,-22),(-24,-23)]:
   sel=(m>=lo)&(m<hi)
   rows.append(dict(branch=branch,published=f.name,Mlo=lo,Mhi=hi,ratio=float(model[sel].sum()/np.power(10.,t['logphi'][sel]).sum())))
  hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
 for z in [.2,.3,.4,.5]:
  for m in [-21.,-22.,-23.]:
   for label,q in [('global',.78),('red',.23),('blue',1.59)]:
    delta=(q-.67)*(z-.1)
    sensitivity.append(dict(branch=branch,z=z,baseline_Mlimit=m,Q=q,population=label,delta_apparent_mag=-delta,conditional_density_ratio=float(cumulative(m+delta)/cumulative(m))))
(O/'RESULTS.json').write_text(json.dumps(dict(comparison=rows,evolution_sensitivity=sensitivity,sha256=hashes,limitations='Reference z=.1, M-5logh and h^3/Mpc^3/mag; bin-integrated cumulative LF. Published total r LF, different K/evolution calibration and survey selection. No covariance/significance claim. Q sensitivities hold LF and colour/K fixed; red/blue Q applied to total LF only to show response scale, not population predictions. Executed wsys/nowsys choice unresolved.'),indent=2)+'\n')
for row in rows:
 if row['branch']=='wsys': print(row)
print('Q sensitivity z=.5 Mlimit=-23:',[r for r in sensitivity if r['branch']=='wsys' and r['z']==.5 and r['baseline_Mlimit']==-23])
