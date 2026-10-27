"""Bounded CPU scaling/GPU resume checks for the authorized survey preflight.

Payload only, no Slurm calls, retry logic, Page writes or production launch.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import loa_survey_runner_20261016 as runner


def replay_root(root, name, cases, include_catalogue=False):
    destination = root/name
    destination.mkdir(exist_ok=False)
    (destination/'INDEX.json').symlink_to(root/'INDEX.json')
    if include_catalogue:
        (destination/'catalogue.h5').symlink_to(root/'catalogue.h5')
    else:
        for case in cases:
            directory = destination/'tiles'/case['id']
            directory.mkdir(parents=True)
            (directory/'prepared').symlink_to(root/'tiles'/case['id']/'prepared')
    return destination


def cpu_scaling(root):
    runner.require_compute()
    index = runner.read(root/'INDEX.json'); cases = index['preflight']
    # Eight four-thread workers fit the allocated 64 logical CPUs. Each prepares
    # two already-tested geometric tiles in an isolated benchmark root.
    destination = replay_root(root, 'cpu_scaling_8', cases, include_catalogue=True)
    started = time.monotonic(); processes=[]
    for worker in range(8):
        batch = cases[worker::8]
        log = (destination/f'worker_{worker}.log').open('x')
        cmd = [sys.executable, '-u', runner.__file__, 'prepare', '--root', str(destination),
               '--tiles', ','.join(case['id'] for case in batch), '--seconds', '1200']
        processes.append((worker, subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), log))
    errors=[]
    for worker, process, log in processes:
        result = process.wait(); log.close()
        if result:
            errors.append({'worker':worker,'returncode':result})
    elapsed = time.monotonic()-started
    if errors:
        runner.atomic_json(destination/'FAILED.json', {'errors':errors,'seconds':elapsed})
        raise RuntimeError('CPU concurrency benchmark failed; see per-tile logs')
    comparisons=[]
    import numpy as np
    for case in cases:
        with np.load(root/'tiles'/case['id']/'prepared/arrays.npz') as original, np.load(destination/'tiles'/case['id']/'prepared/arrays.npz') as replay:
            if set(original.files) != set(replay.files):
                raise ValueError('Prepared replay keys differ')
            for key in original.files:
                np.testing.assert_array_equal(original[key],replay[key])
        comparisons.append(case['id'])
    serial_receipts=[runner.read(p) for p in root.glob('PREPARE_*.json')]
    serial_receipt=next(r for r in serial_receipts if {v['tile'] for v in r['records']}=={c['id'] for c in cases})
    serial_seconds=serial_receipt['seconds']
    record=dict(passed=True,workers=8,tiles=comparisons,seconds=elapsed,
                serial_seconds=serial_seconds,observed_speedup=serial_seconds/elapsed,
                job=os.environ['SLURM_JOB_ID'],includes_worker_startup_and_index_checks=True)
    runner.atomic_json(root/'CPU_SCALING.json',record)
    print(json.dumps(record),flush=True)


def gpu_preflight(root):
    runner.require_compute()
    index = runner.read(root/'INDEX.json'); cases=index['preflight']
    first=cases[0]['id']
    extension=[next(r for r in cases if r['cap']=='NGC' and r['shell']==0 and r['stratum']=='interior'),
               next(r for r in cases if r['cap']=='SGC' and r['shell']==3 and r['stratum']=='edge')]
    started=time.monotonic(); stages=[]

    def run(name, *arguments):
        remaining=6600-(time.monotonic()-started)
        if remaining < 120:
            raise TimeoutError('Preflight work deadline reached; no new allocation or retry')
        command=[sys.executable,'-u',runner.__file__,'infer','--root',str(root),
                 '--seconds',str(remaining-60),*arguments]
        tick=time.monotonic()
        result=subprocess.run(command)
        record=dict(name=name,command=command,seconds=time.monotonic()-tick,returncode=result.returncode)
        stages.append(record)
        runner.atomic_json(root/('STAGE_'+name+'.json'),record)
        if result.returncode:
            raise RuntimeError('Preflight stage failed: '+name)

    run('intentional_stop','--tiles',first,'--max-new','1')
    first_path=root/'tiles'/first/'legacy_v_17/000'
    first_sha=runner.sha(first_path/'arrays.npz')
    first_mtime=(first_path/'arrays.npz').stat().st_mtime_ns
    run('preview')
    if first_sha != runner.sha(first_path/'arrays.npz') or first_mtime != (first_path/'arrays.npz').stat().st_mtime_ns:
        raise ValueError('Resume rewrote a completed draw')
    run('noop_resume','--tiles',first,'--stop','1')
    latest=max(root.glob('INFER_*.json'), key=lambda p:p.stat().st_mtime_ns)
    noop=runner.read(latest)
    if noop['new'] != 0 or noop['reused'] != 6:
        raise ValueError('No-op resume did not reuse all six fits')
    run('extend_32','--tiles',','.join(r['id'] for r in extension),'--first','4','--stop','32')
    replay=replay_root(root,'gpu_replay',[cases[0]])
    run('independent_replay','--root',str(replay),'--tiles',first,'--stop','1')
    import numpy as np
    compared=[]
    for seed in runner.SEEDS:
        paths=[f'coarse_{seed}']+[f'{a}_{p}_{seed}' for a,p in runner.FAMILIES]
        for model in paths:
            relative=Path('tiles')/first/model/'000/arrays.npz'
            with np.load(root/relative) as expected, np.load(replay/relative) as actual:
                if set(expected.files) != set(actual.files):
                    raise ValueError('Replay array keys differ')
                for key in expected.files:
                    np.testing.assert_array_equal(expected[key],actual[key])
            compared.append(str(relative))
    record=dict(passed=True,stages=stages,seconds=time.monotonic()-started,
                preserved_completed_draw=True,noop_resume=dict(new=noop['new'],reused=noop['reused']),
                exact_replay_arrays=compared,extension_tiles=[r['id'] for r in extension],
                job=os.environ['SLURM_JOB_ID'],production_launched=False)
    runner.atomic_json(root/'GPU_PREFLIGHT.json',record)
    summarize(root)


def summarize(root):
    runner.require_compute()
    import numpy as np
    index=runner.read(root/'INDEX.json'); proof=runner.read(root/'GPU_PREFLIGHT.json')
    source=runner.digest(index); records=[]; coarse=[]; tiles=[]
    for case in index['preflight']:
        base=root/'tiles'/case['id']; prep=runner.verified_bundle(base/'prepared',source)
        n=32 if case['id'] in proof['extension_tiles'] else 4
        tile_records=[]; tile_coarse=[]
        for seed in runner.SEEDS:
            for draw in range(n):
                cp=base/f'coarse_{seed}'/f'{draw:03d}'
                cr=runner.read(cp/'COMPLETE.json')
                runner.verified_bundle(cp,cr['binding'])
                if cr['seed'] != seed or cr['draw'] != draw or cr['tile'] != case['id']:
                    raise ValueError('Coarse identity mismatch')
                tile_coarse.append(cr)
                for architecture,prediction in runner.FAMILIES:
                    model=f'{architecture}_{prediction}_{seed}'; path=base/model/f'{draw:03d}'
                    record=runner.verified_bundle(path,runner.digest([cr['binding'],model,draw,cr['sha256']]))
                    if record is None or record['model']!=model or record['draw']!=draw or record['tile']!=case['id']:
                        raise ValueError('Missing/misaddressed draw')
                    with np.load(path/'arrays.npz') as f:
                        for key in f.files:
                            if not np.isfinite(f[key]).all():
                                raise ValueError('Nonfinite output')
                        if f['rho'].shape!=(32,16,16) or f['tensor'].shape!=(32,16,16,6):
                            raise ValueError('Incorrect core geometry')
                        if f['galaxy_eigen'].shape!=(case['retained_rows'],3):
                            raise ValueError('Incorrect galaxy output count')
                        np.testing.assert_allclose(f['tensor'][...,[0,3,5]].sum(-1),f['rho']-1,rtol=3e-6,atol=3e-6)
                        np.testing.assert_allclose(f['galaxy_eigen'].sum(-1),f['galaxy_rho']-1,rtol=3e-6,atol=3e-6)
                        if np.any(f['rho']<=0) or np.any(np.diff(f['galaxy_eigen'],axis=-1)<0):
                            raise ValueError('Invalid saved density/eigenvalue order')
                    tile_records.append(record)
        records.extend(tile_records); coarse.extend(tile_coarse)
        preview_seconds=sum(r['including_tides_seconds'] for r in tile_records if r['draw']<4)+sum(r['seconds'] for r in tile_coarse if r['draw']<4)
        tiles.append(dict(tile=case['id'],cap=case['cap'],shell=case['shell'],stratum=case['stratum'],
                          preview_compute_seconds=preview_seconds,prepare_seconds=prep['seconds'],
                          draws=n,retained_rows=case['retained_rows'],eligible=prep['property_eligible'],
                          prepared_bytes=prep['bytes']))
    strata=[]; estimated_preview=0.; estimated_cpu=0.
    for cap in ('NGC','SGC'):
        for shell in range(4):
            population=sum(r['cap']==cap and r['shell']==shell for r in index['cases'])
            sample=[r for r in tiles if r['cap']==cap and r['shell']==shell]
            mean=float(np.mean([r['preview_compute_seconds'] for r in sample]))
            prep_mean=float(np.mean([r['prepare_seconds'] for r in sample]))
            estimated_preview+=population*mean
            estimated_cpu+=population*prep_mean
            strata.append(dict(cap=cap,shell=shell,population_tiles=population,
                               preview_seconds_per_tile=mean,prepare_seconds_per_tile=prep_mean))
    # Project from the shell-weighted four-draw measurements, then separately
    # report full end-to-end preview wall time (model setup, checks and writes).
    stage=next(s for s in proof['stages'] if s['name']=='preview')
    measured_preview=sum(r['preview_compute_seconds'] for r in tiles)
    # Include the work already committed by the deliberate interrupted stage.
    # Counting its full wall time is conservative: it includes a second startup.
    interrupted=next(s for s in proof['stages'] if s['name']=='intentional_stop')
    overhead_factor=max(1.,(stage['seconds']+interrupted['seconds'])/measured_preview)
    cpu=runner.read(root/'CPU_SCALING.json') if (root/'CPU_SCALING.json').exists() else None
    serial_preps=[runner.read(p) for p in root.glob('PREPARE_*.json')]
    serial_prep=next(r for r in serial_preps if len(r['records'])==16)
    prep_io_factor=serial_prep['seconds']/sum(r['seconds'] for r in serial_prep['records'])
    estimated_cpu*=prep_io_factor
    result=dict(passed=True,fine_draws=len(records),coarse_draws=len(coarse),tiles=tiles,strata=strata,
        preview_compute_gpu_hours=estimated_preview/3600,
        full_32_compute_gpu_hours=estimated_preview*8/3600,
        measured_preview_wall_over_compute=overhead_factor,
        full_32_with_measured_overhead_gpu_hours=estimated_preview*8*overhead_factor/3600,
        serial_prepare_node_hours=estimated_cpu/3600,
        prepare_wall_over_condition_compute=prep_io_factor,
        prepare_node_hours_with_measured_scaling=estimated_cpu/3600/cpu['observed_speedup'] if cpu else None,
        index_seconds=index['seconds'],cpu_scaling=cpu,
        peak_gpu_bytes=max(r['peak_gpu_bytes'] for r in records),
        peak_gpu_reserved_bytes=max(r['peak_gpu_reserved_bytes'] for r in records),
        peak_rss_kib=max(r['peak_rss_kib'] for r in records),
        max_mass_relative=max(r['mass_relative'] for r in records),
        max_core_trace_abs=max(r['trace_maxabs'] for r in records),
        full_parent_retained_rows=sum(r['retained_rows'] for r in index['cases']),
        full_parent_excluded_interpolation_edges=sum(r['excluded_interpolation_edge'] for r in index['cases']),
        closure=runner.CLOSURE,production_authorized=False,
        limitations=['Two geometric tiles per cap/redshift shell, not a precision runtime forecast',
                     'CPU concurrency measured on sixteen tiles, not sustained full-survey I/O',
                     'Final merge/Page packaging budget remains unmeasured',
                     'Numerical validity is not posterior calibration'])
    runner.atomic_json(root/'MEASURED_BUDGET.json',result)
    print(json.dumps(result,indent=2),flush=True)


def budget_closeout(root):
    """Small metadata-only reconciliation; no numerical arrays or GPU work.

    Preserve the first budget receipt and explicitly account for the deliberate
    interruption's work/startup, which was absent from the resumed stage timer.
    """
    import tarfile
    original=runner.read(root/'MEASURED_BUDGET.json')
    proof=runner.read(root/'GPU_PREFLIGHT.json')
    preview=next(r['seconds'] for r in proof['stages'] if r['name']=='preview')
    interrupted=next(r['seconds'] for r in proof['stages'] if r['name']=='intentional_stop')
    factor=original['measured_preview_wall_over_compute']*(preview+interrupted)/preview
    full=original['full_32_compute_gpu_hours']*factor
    result=dict(passed=True,original_budget_sha256=runner.sha(root/'MEASURED_BUDGET.json'),
                preview_stage_seconds=preview,intentional_stop_stage_seconds=interrupted,
                conservative_wall_over_compute=factor,
                preview_with_restart_overhead_gpu_hours=full/8,
                additional_28_with_restart_overhead_gpu_hours=full*7/8,
                full_32_with_restart_overhead_gpu_hours=full,
                core_compute_gpu_hours=original['full_32_compute_gpu_hours'],
                preparation_node_hours=original['prepare_node_hours_with_measured_scaling'],
                preparation_plus_fresh_index_node_hours=original['prepare_node_hours_with_measured_scaling']+original['index_seconds']/3600,
                note='Conservative: includes complete intentional-stop wall time and duplicate startup; not a confidence interval',
                source_sha256=runner.sha(__file__),production_authorized=False)
    archive=root/'budget_closeout_source.tar'
    with tarfile.open(archive,'x') as stream:
        stream.add(__file__,arcname=Path(__file__).name,recursive=False)
    result['source_archive_sha256']=runner.sha(archive)
    runner.atomic_json(root/'BUDGET_CLOSEOUT.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['cpu-scaling','gpu','summarize','budget-closeout'])
    parser.add_argument('--root',type=Path,default=runner.DEFAULT_ROOT)
    args=parser.parse_args()
    {'cpu-scaling':cpu_scaling,'gpu':gpu_preflight,'summarize':summarize,'budget-closeout':budget_closeout}[args.action](args.root)
