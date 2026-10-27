"""Bounded actual-checkpoint numerical parity diagnostic, no posterior products."""
import json,os,time,sys
from pathlib import Path
import numpy as np
import torch
import loa_shared_sixfit_20261026 as s
assert os.environ.get('SLURM_JOB_ID') and torch.cuda.device_count()==1
s.settings.setup(False);torch.backends.cuda.enable_flash_sdp(True);torch.backends.cuda.enable_mem_efficient_sdp(True)
index=s.read(s.INPUT/'INDEX.json');rows=s.select(index,'NGC',False)[:1]
models,hashes,normalizer=s.evaluate.load_models(17,13312,26624);chart=s.views.load_chart(normalizer)
import loa_survey_runner_20261016 as runner
with torch.no_grad():
    obs,_=s.load_observations(rows,chart,runner.digest(index));old,new,wide=obs[0]
    condition=s.Condition(new,wide,old.new_zeros((1,1,48,48,48)),'sampled')
    model=s.FineModel('unet','x').cuda().eval();ck=index['checkpoints']['unet_x_17'];assert s.sha(ck['path'])==ck['sha256']
    state=torch.load(ck['path'],map_location='cpu',weights_only=False);model.load_state_dict(state['ema']);del state
    g=s.geometry([rows[0]['midpoint']],True);noise=s.addressed_noise(g,'NGC',17,'expanded-parity')
    adapter=s.ExpandedAdapter(model);records=[]
    for nfe in (4,128):
        reference=s.sample(model,condition,noise,nfe)
        for name,initial in [('native_arithmetic_consensus',noise)]:
            tested,_=s.integrate_expanded(model,[condition],initial,g,nfe,cache_conditions='--uncached' not in sys.argv)
            diff=(reference-tested).double()
            r=dict(nfe=nfe,variant=name,maxabs=float(diff.abs().max()),rmse=float(diff.square().mean().sqrt()),reference_rms=float(reference.double().square().mean().sqrt()))
            records.append(r);print(json.dumps(r),flush=True)
            assert r['maxabs']<1e-6
    name='PARITY_REPAIRED_UNCACHED.json' if '--uncached' in sys.argv else 'PARITY_REPAIRED.json'
    s.write(Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1')/name,records)
