"""Paired real/RSD and finite-count controls for the three fixed mock regions."""
import argparse,json,os,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter,map_coordinates
from compare_offslab import shape
from mock_truth import ROOT,R
from prepare_parent import digest


def weighted_shape(x,w):
 keep=np.isfinite(w)&(w>0);x=x[keep];w=w[keep];w=w/w.sum();mean=w@x;dx=x-mean
 e=np.maximum(np.linalg.eigvalsh((dx*w[:,None]).T@dx),0)
 return np.array([(e[2]-e[1])/e[2],(e[1]-e[0])/e[2],e[0]/e[2]])


def stat(x):
 return {'n':len(x),'LPS':shape(x)[0].tolist() if len(x)>=3 else None}


def thin(x,rs,n,rng,reps=512):
 if len(x)<n or n<3:return {'available':False},None
 indices=np.stack([rng.choice(len(x),n,replace=False) for _ in range(reps)])
 a=np.array([shape(x[q])[0] for q in indices]);b=np.array([shape(rs[q])[0] for q in indices]);delta=(b[:,1]-b[:,0])-(a[:,1]-a[:,0])
 def summary(v):return {'median_LPS':np.median(v,axis=0).tolist(),'q025_q975_LPS':np.quantile(v,[.025,.975],axis=0).tolist(),'dominant_counts':np.bincount(v.argmax(1),minlength=3).tolist()}
 return {'available':True,'draws':reps,'n':n,'real':summary(a),'rs_same_ids':summary(b),'rs_minus_real_PminusL_quantiles':np.quantile(delta,[.025,.5,.975]).tolist()},(indices,a,b)


def main():
 ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);args=ap.parse_args();out=args.root;assert os.environ.get('SLURM_JOB_ID')
 extraction=json.loads((out/'EXTRACTION.json').read_text())
 if (out/'COMPARISON.json').exists():raise FileExistsError('Completed comparison exists')
 sys.path.insert(0,'/global/u2/d/dkololgi/TNG/Illustris')
 from workflows.abacus_tweb.p12a_coordinate_sample_audit import mapped_member,grid_indices
 from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
 dm=next((ROOT/'ph006/targets/density').glob('*.manifest.json'));meta=json.loads(dm.read_text());density=np.load(meta['build']['output'],mmap_mode='r');assert density.shape==(2048,2048,2048) and meta['build']['processed_file_count']==136
 dx=2000/2048;mean=meta['build']['deposited_count']/2048**3
 tw=ROOT/'ph006/targets/tweb/backend_optimized_ngrid_2048_rsmooth_7';slabs=[]
 for p in sorted(tw.glob('abacus_cactus_tweb_rank*.npz')):
  with np.load(p) as f:slabs.append((int(f['x_start']),int(f['x_end']),p))
 ci=np.load(ROOT/'ph006/p1_canonical/canonical_index.npz');lookup=np.full(int(ci['targetid'].max())+1,-1,dtype='i8');lookup[ci['targetid']]=np.arange(len(ci['targetid']));points=np.load(ROOT/'ph006/p1_canonical/points.npy',mmap_mode='r')
 rng=np.random.default_rng(171019);reports=[]
 for row in extraction['regions']:
  key=row['key'];p=out/(key+'.npz');assert digest(p)==row['sha256'];d=np.load(p);raw=d['raw'];cen=int(d['central_index']);ids=d['targetid'];obs=d['observed'];A=d['transform'];native=d['native_mpc_h'];host=d['host_native']
  assert max(row['transform_residual_Mpc_quantiles'])<.002
  # CutSky coordinates differ from the native snapshot by L/2 in each axis.
  # Validate this fixed offset with 16 nearest valid hosts, independent of labels.
  real=d['real_Mpc']-d['real_Mpc'][cen];rs=d['rs_Mpc']-d['rs_Mpc'][cen]
  candidates=np.flatnonzero((raw['HALO_INDEX']>=0)&(raw['BOX_INDEX']>=0));near=candidates[np.argsort(np.linalg.norm(real[candidates],axis=1))[:16]]
  hostpos=np.empty((len(near),3))
  for file in np.unique(raw['FILE_NUM'][near]):
   path=Path(row['halo_path']).parent/f'halo_info_{int(file):03d}.asdf';cat=CompaSOHaloCatalog(path,fields=['x_com'],subsamples=False,cleaned=False,convert_units=True,verbose=False)
   take=raw['FILE_NUM'][near]==file;hostpos[take]=cat.halos['x_com'][raw['HALO_INDEX'][near[take]]];del cat
  corrected=(native[near]+1000)%2000;offset=(hostpos-corrected+1000)%2000-1000;dist=np.linalg.norm(offset,axis=1)
  # Differences include physical HOD galaxy-host offsets; require the origin
  # correction to reduce the median mismatch from ~sqrt(3)*L/2 to <1 Mpc/h.
  assert np.median(dist)<1
  ih=grid_indices(host[None])[0];covered=0
  for lo,hi,path in slabs:
   if lo<=ih[0]<hi:
    values=mapped_member(path,'eig_vals')[:,ih[0]-lo,ih[1],ih[2]]
    np.testing.assert_allclose(values,row['truth_eigenvalues'],atol=0,rtol=0);covered+=1
  assert covered==1
  origin=(native[cen]+1000)%2000
  base=np.floor(origin/dx).astype(int);offsets=np.arange(-48,49);axes=[(b+offsets)%2048 for b in base]
  patch=np.array(density[np.ix_(*axes)],dtype=float)
  smooth=gaussian_filter(patch,7/dx,mode='nearest',truncate=4)
  line=np.linspace(-1.5*R,1.5*R,97);xx=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1);native_offset=xx@np.linalg.inv(A)
  local=(native_offset+origin-base*dx)/dx+48;mc=map_coordinates(patch,local.reshape(-1,3).T,order=1,mode='nearest').reshape(xx.shape[:3]);sc=map_coordinates(smooth,local.reshape(-1,3).T,order=1,mode='nearest').reshape(xx.shape[:3])
  inside=np.linalg.norm(xx,axis=3)<=R;x=xx[inside];mass=mc[inside];excess=np.maximum(mass-mean,0)
  matter={'total_mass_LPS':weighted_shape(x,mass).tolist(),'positive_excess_LPS':weighted_shape(x,excess).tolist(),'R7_smoothed_mass_LPS':weighted_shape(x,sc[inside]).tolist(),'grid_cell_mpc_h':dx,'native_center':origin.tolist(),'central_galaxy_host_offset_mpc_h':float(np.linalg.norm((host-origin+1000)%2000-1000))}
  # Use the exact observed catalogue positions for observed galaxies.
  orows=lookup[ids[obs]];assert (orows>=0).all();op=points[orows,:3];oriented=np.column_stack((-op[:,1],-op[:,0],op[:,2]));center_obs=np.column_stack((-points[[row['parent_node_id']],1],-points[[row['parent_node_id']],0],points[[row['parent_node_id']],2]))[0]
  exact_rs=np.full_like(rs,np.nan);exact_rs[obs]=oriented-center_obs
  central=np.arange(len(raw))==cen;inreal=(np.linalg.norm(real,axis=1)<=R)&~central;inrs=(np.linalg.norm(exact_rs,axis=1)<=R)&obs&~central
  assert inrs.sum()==row['neighbors_R'];np.testing.assert_allclose(shape(exact_rs[inrs])[0],row['LPS_R'],atol=1e-10,rtol=0)
  observed_real=inreal&obs;bright=inreal&(ids>0);allreal=inreal;n=int(inrs.sum())
  stats={'all_available_galaxies_real':stat(real[allreal]),'bright_parent_real':stat(real[bright]),'observed_fixed_real_ball':stat(real[observed_real]),'observed_fixed_real_ids_in_rs':stat(exact_rs[observed_real]),'observed_fixed_rs_ball':stat(exact_rs[inrs]),'observed_fixed_rs_ids_in_real':stat(real[inrs])}
  controls={};arrays={}
  for name,mask in [('all_available_galaxies',allreal),('bright_parent',bright)]:
   result,saved=thin(real[mask],rs[mask],n,rng);controls[name]=result
   if saved is not None:arrays[name+'_draw_indices'],arrays[name+'_thin_real_LPS'],arrays[name+'_thin_rs_LPS']=saved
  # An unbiased matter sampling control; sampling density is interpolated on the
  # display grid, not a replacement for galaxy bias or the DESI selection function.
  matter_indices=rng.choice(len(x),size=(512,n),replace=True,p=mass/mass.sum());ms=np.array([shape(x[q])[0] for q in matter_indices]);matter['count_matched_sampling']={'n':n,'draws':512,'median_LPS':np.median(ms,axis=0).tolist(),'q025_q975_LPS':np.quantile(ms,[.025,.975],axis=0).tolist(),'dominant_counts':np.bincount(ms.argmax(1),minlength=3).tolist()}
  rec={'key':key,'targetid':row['targetid'],'coordinate_validation':{'heldout_mapping_max_Mpc':row['transform_residual_Mpc_quantiles'][-1],'host_origin_check_rows':len(near),'host_offset_after_origin_correction_quantiles':np.quantile(dist,[0,.5,.95,1]).tolist(),'native_tidal_center_exact_replay':True},'matter':matter,'galaxies':stats,'thinning':controls,'observed_count_target':n,'observed_density_Mpc3':float(n/(4*np.pi*R**3/3)),'rs_real_membership':{'intersection':int((inrs&observed_real).sum()),'rs_only':int((inrs&~observed_real).sum()),'real_only':int((observed_real&~inrs).sum())},'same_rs_ids_delta_PminusL':float((np.array(stats['observed_fixed_rs_ball']['LPS'])-np.array(stats['observed_fixed_rs_ids_in_real']['LPS']))@np.array([-1,1,0]))}
  dest=out/(key+'_comparison.npz');np.savez_compressed(dest,grid_axis_Mpc=line,matter_counts=mc.astype('f4'),matter_R7_counts=sc.astype('f4'),cosmic_mean=mean,real_delta_Mpc=real,rs_delta_Mpc=rs,exact_observed_rs_delta_Mpc=exact_rs,targetid=ids,observed=obs,inreal=inreal,inrs=inrs,central=central,matter_thin_LPS=ms,**arrays)
  rec['arrays_sha256']=digest(dest);reports.append(rec);print(json.dumps(rec),flush=True)
 report={'job':os.environ['SLURM_JOB_ID'],'regions':reports,'seed':171019,'density_manifest_sha256':digest(dm),'density_source':{'path':meta['build']['output'],'bytes':Path(meta['build']['output']).stat().st_size,'mtime_ns':Path(meta['build']['output']).stat().st_mtime_ns,'full_hash_verified':False},'extraction_sha256':digest(out/'EXTRACTION.json'),'script_sha256':digest(__file__),'scope':'Three fixed illustrative regions, already-exposed ph006. Native matter 10% A+B particle TSC grid at 0.977 Mpc/h. Local geometry-only affine transform validated on held-out galaxies; no truth-label fitting. Uniform count-matched thinning is not a DESI selection emulator. Thin real/rs pairs use identical galaxy IDs. Fixed-rs observed cohort replays prior counts/geometry exactly. No model changes.'}
 (out/'COMPARISON.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
