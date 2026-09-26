"""Truth-free Loa versus licensed ph006 input-distribution diagnostics."""
import json,sys
from pathlib import Path
import numpy as np
import healpy as hp
import fitsio
from scipy.stats import ks_2samp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import loa_trial as t

def main(root):
 ready=t.read(root/'INPUTS_READY.json'); d=np.load(root/'catalogue.npz'); z=d['z']; cap=d['cap']; shell=d['shell']
 selection=t.read(ready['selection_path']); galaxy=np.load(root/'galaxy_angular_support.npy'); random=np.load(root/'random_angular.npz'); angles=np.load(root/'boundary_angles.npz')
 pix=np.arange(len(galaxy)); ra,dec=hp.pix2ang(256,pix,lonlat=True); pc=t.galactic_cap(ra,dec)
 mock=t.read(t.B/'ph006/p1_canonical/manifest.json'); mangular=np.load(t.B/'ph006/p3_fields/angular_support_nside256.npz'); print('angular keys',mangular.files,flush=True)
 mcounts=mangular['counts'] if 'counts' in mangular else mangular['count']
 support=np.zeros(len(z),bool); distance=np.zeros(len(z),'f4')
 for c,name in [(0,'SGC'),(1,'NGC')]:
  rows=np.flatnonzero(cap==c)
  distance[rows],support[rows]=t.response_at(d['xyz'][rows],ready['grids'][name],random['support']&((random['domain']//2)==c),angles[name],selection)
 report=[]
 for c,name in [(0,'SGC'),(1,'NGC')]:
  area=(galaxy&(pc==c)).sum()*hp.nside2pixarea(256); ma=((mcounts>0)&(pc==c)).sum()*hp.nside2pixarea(256)
  curve=selection['rotations']['0']['caps'][name]
  for sh,(a,b) in enumerate(zip([.15,.25,.35,.45],[.25,.35,.45,.55])):
   take=(cap==c)&(shell==sh); zg=np.linspace(a,b,1001); rg=np.interp(zg,selection['cosmology']['redshift_grid'],selection['cosmology']['radius_grid_mpc']); ng=np.interp(zg,curve['grid_z'],curve['ntilde']); expected=area*np.trapz(ng*rg**2,rg)
   mn=mock['counts']['by_shell'][f'{a:.2f}_{b:.2f}'][name]
   report.append(dict(cap=name,shell=sh,rows=int(take.sum()),supported=int(support[take].sum()),angular_area_deg2=float(area*(180/np.pi)**2),observed_over_frozen_selection=float(take.sum()/expected),ph006_area_normalized_count_ratio=float(take.sum()/area/(mn/ma)),boundary_quantiles=np.quantile(distance[take],[0,.05,.5,.95,1]).tolist()))
 vac=fitsio.read(root/'DESI_LOA_P12A_HALO48_CANARY_VAC.fits'); x=np.column_stack([vac['BASE_EIGENVALUES'],vac['Z'],np.log(vac['NTILDE_MPC3']),vac['CAP'],np.log1p(vac['BOUNDARY_MPC'])]); ref=np.load(t.C/'dataset/ph006_selection_sample.npz')['context']; comparisons=[]
 for c in [0,1]:
  for a,b in zip([.15,.25,.35,.45],[.25,.35,.45,.55]):
   v=x[(x[:,5]==c)&(x[:,3]>=a)&(x[:,3]<b)&vac['SUPPORTED']]; r=ref[(ref[:,5]==c)&(ref[:,3]>=a)&(ref[:,3]<b)]
   if not len(v) or not len(r):continue
   q=np.quantile(r,[.16,.5,.84],axis=0); scale=np.maximum((q[2]-q[0])/2,1e-8)
   comparisons.append(dict(cap=c,zlo=a,observed_rows=len(v),reference_rows=len(r),median_shift_reference_sigma=((np.median(v,axis=0)-q[1])/scale).tolist(),ks_distance=[float(ks_2samp(v[:,i],r[:,i]).statistic) for i in range(7)]))
 out=dict(schema='loa-input-distributions-v1',selection_comparison=report,canary_feature_comparison=comparisons,science_rows=int((shell>=0).sum()),supported_science_rows=int((support&(shell>=0)).sum()),note='Canary is geometry selected, not representative. KS distances are descriptive, not galaxy-independent significance tests. Area ratios do not model detailed angular response or cosmic variance. No automatic calibration pass or selection refit.',ready_input_sha256=t.digest(root/'INPUTS_READY.json'))
 t.save(root/'INPUT_DISTRIBUTIONS.json',out);print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main(Path(sys.argv[1]))
