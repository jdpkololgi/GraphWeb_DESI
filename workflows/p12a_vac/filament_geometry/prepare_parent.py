"""Fail-closed frozen VAC extraction. Execute full scan on an approved CPU node."""
import argparse, hashlib, json, os
from pathlib import Path
import numpy as np

EXPECTED='cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc'
CENTER=np.array([0.,1025.,0.])
RADIUS=320.

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def coordinates(ra,dec,distance):
    ra=np.deg2rad(ra-180); dec=np.deg2rad(dec)
    return np.column_stack((distance*np.cos(dec)*np.sin(ra),distance*np.cos(dec)*np.cos(ra),distance*np.sin(dec)))

def rotation(angle):
    a=np.deg2rad(angle);c,s=np.cos(a),np.sin(a)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Approved compute allocation required')
    from astropy.io import fits
    from astropy.cosmology import Planck18
    m=json.loads(a.manifest.read_text());source=Path(m['catalogue_path']);rows_path=a.manifest.parent/'plotted_galaxies.npz'
    if m['catalogue_sha256']!=EXPECTED:raise ValueError('Unexpected catalogue binding')
    if a.out.exists():raise FileExistsError('Use an absent isolated output directory')
    if digest(rows_path)!=m['output_hashes']['plotted_galaxies.npz']:raise ValueError('Original rows changed')
    if digest(source)!=EXPECTED:raise ValueError('Frozen catalogue hash mismatch')
    for path,sha in m['model_bindings'].items():
        if digest(path)!=sha:raise ValueError('Model binding changed: '+path)
    original=dict(np.load(rows_path));parts=[];xyz=[];indices=[]
    zgrid=np.linspace(.15,.55,10001);distance=Planck18.comoving_distance(zgrid).value
    columns=['TARGETID','RA','DEC','Z','P_WEB','SUPPORTED','QUALITY','BOUNDARY_MPC','CORE_ID','EIGENVALUE_MEAN','NTILDE_MPC3']
    with fits.open(source,memmap=True) as f:
        for start in range(0,len(f[1].data),200000):
            end=min(start+200000,len(f[1].data));d=f[1].data[start:end]
            z=np.asarray(d['Z']);assert np.all((z>=.15)&(z<.55))
            pos=coordinates(d['RA'],d['DEC'],np.interp(z,zgrid,distance));keep=np.linalg.norm(pos-CENTER,axis=1)<=RADIUS
            parts.append({k:np.array(d[k][keep]) for k in columns});xyz.append(pos[keep]);indices.append(np.arange(start,end)[keep])
    data={k:np.concatenate([b[k] for b in parts]) for k in columns};data['xyz_Mpc']=np.concatenate(xyz);data['catalogue_row']=np.concatenate(indices)
    assert len(np.unique(data['TARGETID']))==len(data['TARGETID'])
    ii=np.searchsorted(data['catalogue_row'],original['catalogue_row'])
    assert np.all(ii<len(data['TARGETID']))
    for k in ['TARGETID','RA','DEC','Z','QUALITY','SUPPORTED']:np.testing.assert_array_equal(data[k][ii],original[k])
    np.testing.assert_array_equal(data['P_WEB'][ii,2],original['P_filament'])
    for j,k in enumerate(['x_Mpc','y_Mpc','height_Mpc']):np.testing.assert_allclose(data['xyz_Mpc'][ii,j],original[k],atol=1e-8,rtol=0)
    # All parent tracers retained regardless of posterior support, quality or P_WEB.
    yes=data['SUPPORTED'];assert np.all(np.isfinite(data['P_WEB'][yes]));np.testing.assert_allclose(data['P_WEB'][yes].sum(1),1,atol=1e-6)
    assert digest(source)==EXPECTED
    a.out.mkdir(parents=True);data['original_target_parent_index']=ii
    np.savez_compressed(a.out/'parent.npz',**data)
    config={'center_Mpc':CENTER.tolist(),'parent_radius_Mpc':RADIUS,'full_depths_Mpc':[5,10,20,40,80], 'rotation_about_radial_y_degrees':[0,30,60,90],'normal_offsets_Mpc':[-5,0,5], 'half_widths_Mpc':[180,225], 'radii_Mpc':[.5*7/.6766,7/.6766,1.5*7/.6766], 'min_neighbors':[6,8,12], 'shape_cut_frozen':.5847730874520835,'low_Pfil_cut':.2,'note':'Fixed original target IDs tracked; targets outside each newly selected slab are excluded and counted, never silently replaced. Geometry/count eligibility must be intersected for paired tests. Parent-padding completeness is not survey-mask completeness. VAC parent is not automatically identical to full inference context tracer sample.'}
    assert np.sqrt(180**2+225**2+45**2)+max(config['radii_Mpc'])<RADIUS
    receipt={'stage':'parent extraction only; no scientific verdict','catalogue_sha256':EXPECTED,'manifest_sha256':digest(a.manifest),'script_sha256':digest(__file__),'model_bindings_verified':m['model_bindings'],'parent_sha256':digest(a.out/'parent.npz'),'rows':len(data['TARGETID']),'original_targets':len(ii),'job':os.environ['SLURM_JOB_ID'],'config':config}
    (a.out/'READY.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(receipt))
if __name__=='__main__':main()
