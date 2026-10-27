"""Audit the complete four-draw preview; never authorize extension or submit jobs.

Numerical checks are implementation requirements. Field summaries are descriptive,
equal-tile diagnostics, not calibration tests or independent-galaxy error bars.
Run on an approved CPU node, after every preview shard has completed.
"""
from concurrent.futures import ThreadPoolExecutor
import io
import os
import tarfile
import time

import numpy as np
import loa_survey_production_20261016 as p
r = p.r


def check_arrays(arrays, rows, coarse_core=None):
    shapes = dict(rho=(32,16,16), tensor=(32,16,16,6),
                  galaxy_eigen=(rows,3), galaxy_rho=(rows,))
    if set(arrays) != set(shapes):
        raise ValueError('Unexpected saved field members')
    for key, shape in shapes.items():
        if arrays[key].shape != shape or not np.isfinite(arrays[key]).all():
            raise ValueError('Wrong shape/nonfinite field: '+key)
    if np.any(arrays['rho'] <= 0) or np.any(arrays['galaxy_rho'] <= 0):
        raise ValueError('Nonpositive density')
    if np.any(np.diff(arrays['galaxy_eigen'],axis=-1) < 0):
        raise ValueError('Unordered eigenvalues')
    np.testing.assert_allclose(arrays['tensor'][...,[0,3,5]].sum(-1),
                               arrays['rho']-1,rtol=2e-5,atol=2e-5)
    np.testing.assert_allclose(arrays['galaxy_eigen'].sum(-1),
                               arrays['galaxy_rho']-1,rtol=2e-5,atol=2e-5)
    if coarse_core is not None:
        mean=arrays['rho'].astype(float).reshape(8,4,4,4,4,4).mean(axis=(1,3,5))
        np.testing.assert_allclose(mean,coarse_core,rtol=2e-6,atol=1e-8)


def normalizer_scale(root):
    import torch
    from workflows.sbi import e2e_coupled_views as views
    index=r.read(root/'INDEX.json'); hashes=[]
    for seed in r.SEEDS:
        record=index['checkpoints'][f'coarse_{seed}']
        if r.sha(record['path'])!=record['sha256']:
            raise ValueError('Coarse checkpoint changed before normalization audit')
        state=torch.load(record['path'],map_location='cpu',weights_only=False)
        hashes.append(state['binding']['normalizer'])
        del state
    if hashes[0]!=hashes[1]:
        raise ValueError('Expected common frozen normalizer')
    value=views.load_chart(hashes[0])['coarse_logrho']
    return value['mean'][0],value['std'][0]


def audit_tile(root, row, scale):
    tile=row['id']; meta=p.verify_archive(root,tile,4)
    with tarfile.open(p.archive_path(root,tile,4)/'tile.tar','r') as archive:
        def load(name, keys=None):
            with archive.extractfile(name) as stream, np.load(io.BytesIO(stream.read())) as saved:
                return {key:saved[key] for key in (saved.files if keys is None else keys)}
        keys=('TARGETID','Z','SUPPORTED','core_support','CIGALE_MATCHED','CIGALE_AMBIGUOUS',
              'CIGALE_SPECTYPE','MASS_CG_15','r_bright')
        prep=load('prepared/arrays.npz',keys)
        ids=prep['TARGETID']; n=len(ids)
        if n != row['retained_rows'] or len(np.unique(ids)) != n:
            raise ValueError('Prepared ownership/count mismatch')
        mass=prep['MASS_CG_15']
        eligible=(prep['SUPPORTED'].astype(bool)&(prep['core_support']>=.5)&
                  prep['CIGALE_MATCHED'].astype(bool)&~prep['CIGALE_AMBIGUOUS'].astype(bool)&
                  (np.char.strip(prep['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')&
                  np.isfinite(mass)&(mass>=1e9)&(mass<1e12)&(prep['r_bright']<19.5))
        with archive.extractfile('prepared/COMPLETE.json') as stream:
            import json
            receipt=json.load(stream)
        if int(eligible.sum()) != receipt['property_eligible']:
            raise ValueError('Eligibility reconciliation failed')
        coarse_cores={}
        for seed in r.SEEDS:
            for draw in range(4):
                z=load(f'coarse_{seed}/{draw:03d}/arrays.npz')['z']
                if z.shape!=(1,1,48,48,48) or not np.isfinite(z).all():
                    raise ValueError('Invalid stored coarse draw')
                wide=np.exp(z[0,0].astype(float)*scale[1]+scale[0])
                if not np.isfinite(wide).all() or np.any(wide<=0):
                    raise ValueError('Invalid decoded coarse density')
                coarse_cores[seed,draw]=wide[20:28,22:26,22:26]
        summaries={}
        for seed in r.SEEDS:
            for architecture,prediction in r.FAMILIES:
                model=f'{architecture}_{prediction}_{seed}'
                values=[]; classes=np.zeros(4,dtype=np.int64); selected=np.zeros(4,dtype=np.int64)
                for draw in range(4):
                    arrays=load(f'{model}/{draw:03d}/arrays.npz')
                    check_arrays(arrays,n,coarse_cores[seed,draw])
                    rho=arrays['rho'].astype(float)
                    values.append([rho.mean(),rho.std(),rho.min(),rho.max()])
                    labels=(arrays['galaxy_eigen']>.2).sum(-1)
                    classes+=np.bincount(labels,minlength=4)
                    selected+=np.bincount(labels[eligible],minlength=4)
                summaries[model]=dict(draw_core_stats=values,retained_class_draw_counts=classes.tolist(),
                                      eligible_class_draw_counts=selected.tolist())
    return dict(tile=tile,cap=row['cap'],tile_median_z_shell=row['shell'],retained=n,
        eligible=int(eligible.sum()),retained_by_galaxy_z_shell=np.histogram(prep['Z'],r.SHELLS)[0].tolist(),
        eligible_by_galaxy_z_shell=np.histogram(prep['Z'][eligible],r.SHELLS)[0].tolist(),
        archive_bytes=meta['bytes'],archive_sha256=meta['archive_sha256'],models=summaries),ids


def main():
    r.require_compute(); root=p.ROOT; p.check_campaign(root)
    validated=r.read(root/'PREVIEW_AUDITOR_VALIDATED.json')
    if not validated['passed'] or any(r.sha(path)!=value for path,value in validated['sources'].items()):
        raise ValueError('Preview auditor changed since representative validation')
    start=time.monotonic(); index=r.read(root/'INDEX.json')
    expected=[row['id'] for row in index['cases']]
    if len(expected)!=2472 or len(set(expected))!=2472:
        raise ValueError('Unexpected tile ownership')
    found=[]
    for shard in range(103):
        receipt=r.read(root/'shards'/f'{shard:03d}_04.json')
        if (not receipt['passed'] or receipt['draws']!=4 or
                receipt['campaign_sha256']!=r.sha(root/'CAMPAIGN.json')):
            raise ValueError('Bad shard receipt')
        found.extend(receipt['tiles'])
    if found!=expected:
        raise ValueError('Missing/duplicate/out-of-order shards')
    scale=normalizer_scale(root)
    with ThreadPoolExecutor(max_workers=8) as pool:
        results=list(pool.map(lambda row:audit_tile(root,row,scale),index['cases']))
    tiles=[result[0] for result in results]
    ids=np.concatenate([result[1] for result in results])
    if len(ids)!=4628192 or len(np.unique(ids))!=len(ids):
        raise ValueError('Global TARGETID ownership mismatch')
    models={}
    for model in tiles[0]['models']:
        values=np.asarray([x['models'][model]['draw_core_stats'] for x in tiles])
        models[model]=dict(
            core_stats=['mean_rho','sd_rho','minimum_rho','maximum_rho'],
            quantiles=[0,.01,.25,.5,.75,.99,1],
            equal_tile_draw_quantiles=np.quantile(values.reshape(-1,4),[0,.01,.25,.5,.75,.99,1],axis=0).tolist(),
            retained_class_draw_counts=np.sum([x['models'][model]['retained_class_draw_counts'] for x in tiles],axis=0).tolist(),
            eligible_class_draw_counts=np.sum([x['models'][model]['eligible_class_draw_counts'] for x in tiles],axis=0).tolist())
    result=dict(passed=True,extension_authorized=False,requires_scientific_and_budget_review=True,
        campaign_sha256=r.sha(root/'CAMPAIGN.json'),audit_source_sha256=r.sha(__file__),
        tiles=len(tiles),fine_draws=2472*6*4,shared_coarse_draws=2472*2*4,
        retained=len(ids),property_eligible=sum(x['eligible'] for x in tiles),
        retained_by_galaxy_z_shell=np.sum([x['retained_by_galaxy_z_shell'] for x in tiles],axis=0).tolist(),
        eligible_by_galaxy_z_shell=np.sum([x['eligible_by_galaxy_z_shell'] for x in tiles],axis=0).tolist(),
        archive_bytes=sum(x['archive_bytes'] for x in tiles),empty_retained_tiles=sum(x['retained']==0 for x in tiles),
        models=models,per_tile=tiles,seconds=time.monotonic()-start,job=os.environ['SLURM_JOB_ID'],
        interpretation='four-draw numerical/coverage audit only; no truth, calibration or independent-survey posterior claim')
    p.check_campaign(root)
    r.atomic_json(root/'PREVIEW_AUDIT.json',result)
    print({k:v for k,v in result.items() if k not in ('per_tile','models')},flush=True)


if __name__=='__main__': main()
