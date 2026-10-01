"""Self-contained local 3D Loa viewer; display thinning never alters analysis.

Full atlas is a mosaic of local conditional marginals, NOT a joint posterior.
No remote scripts, uploads, fabricated positions, or interpolated tile bridges.
"""
import argparse
import base64
import hashlib
import json
import os
import sys
from pathlib import Path
import numpy as np
from plotly.offline import get_plotlyjs


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def pack(a):
    return base64.b64encode(np.asarray(a,dtype='<f4').tobytes()).decode('ascii')


def sample_rows(ids,maximum):
    """Stable integer hash, with deterministic exact top-k if necessary."""
    x=np.asarray(ids,dtype=np.uint64).copy()
    x^=x>>np.uint64(30);x*=np.uint64(0xbf58476d1ce4e5b9)
    x^=x>>np.uint64(27);x*=np.uint64(0x94d049bb133111eb);x^=x>>np.uint64(31)
    return np.sort(np.argpartition(x,maximum-1)[:maximum]) if len(x)>maximum else np.arange(len(x))


def pooled(a):
    shape=a.shape
    return a.reshape(*shape[:-3],shape[-3]//2,2,shape[-2]//2,2,shape[-1]//2,2).mean(axis=(-1,-3,-5))


def xyz_payload(x):
    return {axis:pack(x[:,i]) for i,axis in enumerate('xyz')}


def build(a):
    if 'SLURM_JOB_ID' not in os.environ:raise RuntimeError('Build on compute node')
    a.out.mkdir(parents=True,exist_ok=False)
    if a.surface_library:sys.path.insert(0,str(a.surface_library))
    import skimage
    from skimage.measure import marching_cubes
    plan=json.loads((a.root/'PLAN.json').read_text())
    if digest(a.root/'vac_properties.npz')!=plan['properties_sha256']:raise ValueError('property hash drift')
    with np.load(a.root/'vac_properties.npz') as f:
        xyz=f['xyz'];caps=f['CAP'];ids=f['TARGETID'];z=f['Z']
    payload={};audit={'plan_sha256':digest(a.root/'PLAN.json'),'views':{},'global_joint_posterior':False,
        'script_sha256':digest(__file__),'template_sha256':digest(Path(__file__).with_name('loa_field_viewer.html')),
        'surface_library_version':skimage.__version__,'display_only':True,'remote_dependencies':False}
    for cap,capid in [('NGC',1),('SGC',0)]:
        rows=np.flatnonzero(caps==capid);picked=rows[sample_rows(ids[rows],100000)]
        overview=dict(**xyz_payload(xyz[picked]),zred=pack(z[picked]),count=len(picked),total=len(rows))
        payload[cap+' overview']=dict(kind='overview',gal=overview,
            note=f'{len(picked):,} of {len(rows):,} VAC galaxies; deterministic TARGETID-hash thinning. Actual equatorial comoving coordinates.')
        receipt=json.loads((a.root/'details'/('detail_'+cap)/'COMPLETE.json').read_text())
        for item in receipt['outputs']:
            if digest(item['path'])!=item['sha256']:raise ValueError('detail hash drift')
        with np.load(receipt['outputs'][0]['path']) as f:r=pooled(f['rho']);support=pooled(f['support'])
        with np.load(receipt['outputs'][1]['path']) as f:gal=f['xyz']
        ngal=len(gal);gal=gal[sample_rows(np.arange(ngal),30000)]
        grid=plan['grids'][cap];mid=np.asarray(receipt['case']['midpoint'])
        origin=np.asarray(grid['origin_mpc_h'])+(mid-[64,48,48])*3.383
        pos=np.stack(np.meshgrid(*[origin[i]+(np.arange(n)+.5)*13.532 for i,n in enumerate(r.shape[1:])],indexing='ij'),-1).reshape(-1,3)
        # Mask unsupported display voxels in the scalar, never interpolate across them.
        good=support.ravel()>=.5
        values=np.log10(np.maximum(r.reshape(len(r),-1),1e-8));values[:,~good]=np.nan
        mean=np.log10(np.maximum(r.mean(0).ravel(),1e-8));mean[~good]=np.nan
        meshes=[]
        for volume in np.concatenate([r.mean(0,keepdims=True),r],axis=0):
            surfaces=[]
            for level in (1.,2.,4.):
                if not volume.min()<level<volume.max():surfaces.append(None);continue
                try:
                    verts,faces,_,_=marching_cubes(volume,level=level,spacing=(13.532,)*3,
                        mask=(support>=.5),allow_degenerate=False)
                except RuntimeError:surfaces.append(None);continue
                verts+=origin+6.766
                surfaces.append(dict(**xyz_payload(verts),i=pack(faces[:,0]),j=pack(faces[:,1]),k=pack(faces[:,2]),level=level))
            meshes.append(surfaces)
        payload[cap+' detail']=dict(kind='detail',**xyz_payload(pos),values=[pack(mean)]+[pack(v) for v in values],
            meshes=meshes,
            gal=dict(**xyz_payload(gal),count=len(gal),total=ngal),
            note=f'16 local conditional draws; 13.532 Mpc/h display cells (2× pooling only). Density shown where pooled support ≥ 0.5. {len(gal):,} of {ngal:,} observed galaxies shown by deterministic row-hash thinning; survey selection remains nonuniform.')
        audit['views'][cap+' detail']=dict(receipt_sha256=digest(a.root/'details'/('detail_'+cap)/'COMPLETE.json'),
            display_cells=len(pos),supported_display_cells=int(good.sum()),galaxies=len(gal),total_galaxies=ngal,draws=len(r),
            legacy_density_replay=receipt['legacy_density_replay'])
    if not a.details_only:
        atlas=json.loads((a.root/'ATLAS_COMPLETE.json').read_text())
        if not atlas['passed'] or atlas['tiles']!=len(plan['tiles']):raise ValueError('atlas incomplete')
        for cap in ('NGC','SGC'):
            positions=[];values=[];tileids=[];tilemeta=[]
            for tid,case in enumerate(plan['tiles']):
                if case['cap']!=cap:continue
                receipt=json.loads((a.root/'tiles'/case['id']/'COMPLETE.json').read_text())
                item=receipt['outputs'][0]
                if digest(item['path'])!=item['sha256']:raise ValueError('atlas tile hash drift')
                with np.load(item['path']) as f:r=pooled(f['rho']);support=pooled(f['support'])
                # One reproducible spatial sublattice of the pooled owned core.
                index=np.indices(support.shape);keep=(support>=.5)&((index.sum(0)%4)==0)
                ii=np.stack([x[keep] for x in index],1)
                origin=np.asarray(plan['grids'][cap]['origin_mpc_h'])+(np.asarray(case['midpoint'])-[32,16,16])*3.383
                positions.append(origin+(ii+.5)*13.532)
                vv=np.concatenate([r.mean(0,keepdims=True),r],axis=0)[:,keep]
                values.append(np.log10(np.maximum(vv,1e-8)));tileids.extend([tid]*len(ii));tilemeta.append(case['id'])
            pos=np.concatenate(positions);vv=np.concatenate(values,axis=1)
            payload[cap+' atlas']=dict(kind='atlas',**xyz_payload(pos),values=[pack(v) for v in vv],tile=pack(tileids),
                gal=payload[cap+' overview']['gal'],
                note=f'{len(tilemeta):,} independently sampled owned tiles; {len(pos):,} displayed density voxels. 13.532 Mpc/h display pooling + 1-in-4 spatial thinning; support ≥ 0.5. Same-index samples across tiles are NOT a global joint draw. No smoothing across tile seams.')
            audit['views'][cap+' atlas']=dict(tiles=len(tilemeta),display_voxels=len(pos),draws=8)
    template=Path(__file__).with_name('loa_field_viewer.html').read_text()
    html=template.replace('__PLOTLY__',get_plotlyjs()).replace('__PAYLOAD__',json.dumps(payload,separators=(',',':')))
    (a.out/'Loa_3D.html').write_text(html)
    audit['html_sha256']=digest(a.out/'Loa_3D.html');audit['bytes']=(a.out/'Loa_3D.html').stat().st_size
    (a.out/'VIEWER.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--surface-library',type=Path)
    p.add_argument('--details-only',action='store_true');build(p.parse_args())
