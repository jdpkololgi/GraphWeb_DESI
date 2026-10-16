"""Independent filter integration, source parity and result checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
import h5py
from scipy.interpolate import interp1d
from astropy.cosmology import Planck18
sys.path.insert(0,'/global/common/software/desi/perlmutter/desiconda/current/code/speclite/v1.0.0/lib/python3.13/site-packages')
from speclite import filters
R=Path('/global/u2/d/dkololgi/fmsu_tokyo_20261010');O=R/'tng_bright_refined_20261015'
p=np.load(O/'inputs.npz');o=np.load(O/'observed_photometry.npz');r=json.loads((O/'RESULTS.json').read_text())
knots=np.array([f.effective_wavelength.value for f in filters.load_filters('sdss2010-g','sdss2010-r','sdss2010-i','sdss2010-z')]);fs=filters.load_filters('decam2014-g','decam2014-r','decam2014-z')
ii=np.linspace(0,len(p['ids'])-1,32,dtype=int);wave=np.arange(2900.,11501.,5.);errs=[]
for j in [0,4,9]:
 z=.205+.01*j
 logs=interp1d(np.log(knots),o['logfnu_nodes'][ii],axis=1,fill_value='extrapolate')(np.log(wave/(1+z)))
 # Direct f_lambda AB integration using speclite's implementation, not our f_nu weights.
 flam=10**logs*3.631e-20*2.99792458e18/wave**2*(1+z)
 mag=fs.get_ab_magnitudes(flam,wave)
 vals=np.array([mag[n] for n in mag.colnames]).T+Planck18.distmod(z).value
 errs.extend(np.abs(vals-o['mag_grz'][j,ii]).ravel())
assert max(errs)<.005
# Fresh source parity on a nonempty sample, including cumulative Subfind offsets.
offset=0;source_rows=0
for part in range(600):
 with h5py.File(f'/pscratch/sd/d/dkololgi/TNG300-1/groups_099/fof_subhalo_tab_099.{part}.hdf5') as f:
  n=int(f['Header'].attrs['Nsubgroups_ThisFile']);ix=np.flatnonzero((p['subhalo_id']>=offset)&(p['subhalo_id']<offset+n));sid=p['subhalo_id'][ix]-offset
  if len(ix):
   assert np.array_equal(f['Subhalo/SubhaloStellarPhotometrics'][:,4:8][sid],p['mags'][ix])
   assert np.allclose(f['Subhalo/SubhaloHalfmassRadType'][:,4][sid]/float(f['Header'].attrs['HubbleParam']),p['rhalf_pkpc'][ix])
   source_rows+=len(ix)
  offset+=n
 if source_rows>=256:break
assert source_rows>=256
old=json.loads((R/'tng_bright_selection_20261015/RESULTS.json').read_text())['results']['high']['curves']
for new,prev in [('uncut','uncut'),('old_powerlaw','bright_color'),('loa','loa')]:
 assert np.allclose(r['results']['0.25_high']['curves'][new]['mean'],old[prev]['mean'],atol=1e-12)
for name,v in r['results'].items():
 assert abs(sum(v['weights'])-1)<1e-12
 for arm,c in v['curves'].items():
  assert np.all((np.array(c['mean'])>=0)&(np.array(c['mean'])<=1))
  if arm!='loa':assert np.isclose(c['mean_gap_pp'],100*np.mean(np.array(c['mean'])-v['curves']['loa']['mean']))
filterroot=Path('/global/common/software/desi/perlmutter/desiconda/current/code/speclite/v1.0.0/lib/python3.13/site-packages/speclite/data/filters')
files=[filterroot/(n+'.ecsv') for n in ['sdss2010-g','sdss2010-r','sdss2010-i','sdss2010-z','decam2014-g','decam2014-r','decam2014-z']]
files+=[Path('/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/lib/python3.11/site-packages/desitarget/cuts.py')]
v={'independent_speclite_5A_integrations':len(errs),'maximum_mag_difference':float(max(errs)),'fresh_source_rows':source_rows,'previous_high_mass_curve_parity':3,'cell_mean_checks':r['checks'],'minimum_valid_spatial_bootstraps':min(c['valid_bootstrap_min'] for v in r['results'].values() for c in v['curves'].values()),'filter_and_selection_source_hashes':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}
(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v,indent=2))
