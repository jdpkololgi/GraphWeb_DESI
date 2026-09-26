"""Paired truth-known density perturbations on eight previously exposed cores."""
import sys,json,copy,dataclasses,os
from pathlib import Path
import numpy as np,torch
import loa_trial as t
from workflows.abacus_tweb.p10_training_contract import P10PhaseBalancedLoader
from workflows.abacus_tweb import p8_train_unet_patch as u
from workflows.abacus_tweb.p8_deterministic_common import increments_to_eigenvalues,unscale_increments
from workflows.sbi.p12a_blind_inference import reconstruct_fmpe
from workflows.sbi.p12_train_base_response_fmpe import sample_posterior,theta_to_eigenvalues
OUT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_selection_sensitivity_20260925_v1')
def predict(model,patch,ck):
 x,p=u.model_inputs(patch,ck['normalization'],'cuda')
 with torch.inference_mode():v=model.head(model.sample_latent(x,p)).cpu().numpy()
 return increments_to_eigenvalues(unscale_increments(v,ck['scaler'])).astype('f4')
def main():
 OUT.mkdir(exist_ok=False);torch.set_num_threads(8);torch.backends.cudnn.allow_tf32=True;torch.backends.cuda.matmul.allow_tf32=False
 sr=t.read(t.C/'summaries/ph006/OOF_SUMMARY_COMPLETE.json');loader=P10PhaseBalancedLoader(Path(sr['contract_root']),include_blind=False);adapter=loader.field_adapter('ph006');p3=t.read(adapter.manifest['p3_manifest']);selection=t.read(Path(sr['contract_root'])/'transforms/field/selection_manifest.json');half=copy.deepcopy(selection)
 for v in half['rotations']['0']['caps'].values():v['ntilde']=(np.array(v['ntilde'])*.5).tolist()
 idx=np.load(p3['canonical_index']);points=np.load(p3['points'],mmap_mode='r');parent=np.load(sr['arrays']['parent_node_id']);order=np.argsort(parent);old=np.load(sr['arrays']['base_prediction'],mmap_mode='r');truth=np.load(sr['arrays']['truth'],mmap_mode='r');resp=np.load(sr['arrays']['response']);cache=t.B/'p12a_random_support_parent_cache_v2/ph006';dist=np.load(cache/'distance_to_support_boundary_mpc.npy');support=np.load(cache/'support_random.npy');schema=t.read(p3['frozen_schema']);spline=t.read(p3['ntilde_spline']);angular=np.load(p3['angular_support']['path'])['support']
 ck=torch.load(sr['checkpoint'],map_location='cuda',weights_only=False);model=u.UPatch().cuda().eval();model.load_state_dict(ck['state_dict']);posterior,pck=reconstruct_fmpe(t.C/'posterior/fmpe_estimator.pt','cuda')
 rng=np.random.default_rng(20260925);keep=rng.random(len(points))<.5;chosen=t.read(t.C/'GOLDEN_OBSERVER_REPLAY_20260925_v3.json')['runs'];runs=[];records=[]
 for meta in chosen:
  core=meta['core'];cap=meta['cap'];name='NGC' if cap else 'SGC';patch=adapter.extract(core,48,u.CHANNELS,alignment_voxels=8);grid=adapter.manifest['caps'][name];frac=(points[:,:3]-grid['origin_mpc'])/grid['cell_mpc']-.5;rows=np.flatnonzero((idx['cap']==cap)&idx['context']&np.all((frac>=patch.context_start-1)&(frac<patch.context_stop),axis=1));at=order[np.searchsorted(parent[order],patch.authoritative_parent_id)];assert np.array_equal(parent[at],patch.authoritative_parent_id)
  eligible=np.flatnonzero(keep[parent[at]]&support[at]);sample=np.sort(rng.choice(eligible,min(64,len(eligible)),replace=False));atq=at[sample];tr=np.array(truth[atq]);ref=predict(model,patch,ck);assert np.allclose(ref,old[at],atol=1e-5,rtol=1e-5)
  baseline_q=None;baseline_p=None
  for arm,thin,sel,factor in [('baseline',False,selection,1.),('half_expected',False,half,.5),('thin_wrong_expected',True,selection,1.),('thin_matched_expected',True,half,.5)]:
   rr=rows[keep[rows]] if thin else rows;fields=t.obs.rebuild_fields(points[rr,:3],idx['shell'][rr],grid,patch.context_start,patch.context_stop,angular,p3['angular_support']['nside'],schema,spline,sel,0,name);values=np.stack([fields[n] for n in patch.channel_names]);new=dataclasses.replace(patch,values=values);base=predict(model,new,ck)[sample]
   if arm=='baseline':
    assert np.allclose(values,patch.values,atol=1e-5,rtol=1e-5);assert np.allclose(base,ref[sample],atol=1e-5,rtol=1e-5)
   x=np.column_stack([base,resp['redshift'][atq],np.log(resp['ntilde_mpc3'][atq]*factor),resp['cap'][atq],np.log1p(dist[atq])]);x=((x-pck['context_mean'])/pck['context_std']).astype('f4');torch.manual_seed(20260925+core);torch.cuda.manual_seed_all(20260925+core);ss=sample_posterior(posterior,x,512,64,'cuda');eig=theta_to_eigenvalues(ss*pck['theta_std']+pck['theta_mean']);q=np.quantile(eig,[.05,.16,.5,.84,.95],axis=1).transpose(1,0,2);cl=(eig>.2).sum(-1);p=np.stack([(cl==j).mean(1) for j in range(4)],1);tc=(tr>.2).sum(1)
   if arm=='baseline':baseline_q=q;baseline_p=p
   out=dict(core=core,cap=cap,core_shell=meta['shell'],arm=arm,rows=len(atq),coverage68=((tr>=q[:,1])&(tr<=q[:,3])).mean(0).tolist(),coverage90=((tr>=q[:,0])&(tr<=q[:,4])).mean(0).tolist(),mean_prob=p.mean(0).tolist(),truth_fraction=np.bincount(tc,minlength=4).astype(float).__truediv__(len(tc)).tolist(),brier=float(((p-np.eye(4)[tc])**2).sum(1).mean()),mean_prob_shift=(p-baseline_p).mean(0).tolist(),median_shift_over_baseline_width=np.median((q[:,2]-baseline_q[:,2])/np.maximum(baseline_q[:,3]-baseline_q[:,1],1e-6),axis=0).tolist());runs.append(out)
   np.savez(OUT/f'core{core}_{arm}.npz',parent_id=parent[atq],redshift=resp['redshift'][atq],truth=tr,base=base,quantiles=q,p_web=p);print(json.dumps(out),flush=True)
 loader.close();t.save(OUT/'RESULTS.json',dict(runs=runs,job=os.environ.get('SLURM_JOB_ID'),script_sha256=t.digest(__file__),baseline_replay_pass=True,protocol=str(Path(__file__).resolve().parents[2]/'docs/p12a_selection_test_protocol_20260925.md'),note='Paired geometry-selected development diagnostic; fixed angular response, not global calibration or a validated DESI correction.'))
if __name__=='__main__':main()
