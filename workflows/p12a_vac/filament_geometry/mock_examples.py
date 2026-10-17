"""Three fixed-rule mock examples; full neighbour truth and frozen FMPE predictions.

Central probabilities are the original saved draws. Neighbour probabilities are
new draws from the same frozen model and saved encoder summaries, never refits.
"""
import argparse,json,os,sys
from pathlib import Path
import numpy as np
from prepare_parent import digest
from compare_offslab import shape,contrast2
from mock_truth import ROOT,CAND,R,CUT,eigen


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Approved compute allocation required')
 if args.out.exists():raise FileExistsError(args.out)
 import torch
 torch.set_num_threads(8);torch.manual_seed(171018)
 evidence=Path(__file__).resolve().parents[3]/'docs/evidence/filament_mock_truth_20261017'
 report=json.loads((evidence/'RESULTS.json').read_text());d=np.load(evidence/'geometry.npz')
 assert digest(evidence/'geometry.npz')==report['geometry_sha256']
 t=(d['truth_eigenvalues']>.2).sum(1);p=d['probability'];g=d['local_1_full']
 cohort=(d['local_projected']>=CUT)&(p[:,2]<.2)&(p.argmax(1)==1)&d['local_1_eligible']&(d['local_1_nfull']>=8)&np.isfinite(g).all(1)
 specs=[('true_wall_line','True wall with a string-like tracer neighbourhood',cohort&(t==1)&(g.argmax(1)==0)),('true_wall_sheet','True wall with a sheet-dominated tracer neighbourhood',cohort&(t==1)&(g.argmax(1)==1)),('filament_error','True filament assigned highest probability to wall',cohort&(t==2))]
 chosen=[]
 for key,title,mask in specs:
  rows=np.flatnonzero(mask);rows=rows[np.argsort(d['parent_node_id'][rows])]
  assert len(rows)>0
  chosen.append((key,title,int(rows[len(rows)//2]),len(rows)))
 phase=ROOT/'ph006/p1_canonical';m=json.loads((phase/'manifest.json').read_text())
 points_path=Path(m['points']);index_path=Path(m['index'])
 assert digest(points_path)==report['sha256'][str(points_path)] and digest(index_path)==report['sha256'][str(index_path)]
 points=np.load(points_path,mmap_mode='r');ci=np.load(index_path);truth_path=ROOT/'training_contract/phases/ph006/parent_eigenvalues.npy'
 assert digest(truth_path)==report['sha256'][str(truth_path)]
 native=np.load(truth_path,mmap_mode='r');from astropy.cosmology import Planck18
 bounds=Planck18.comoving_distance([.15,.55]).value
 rad=np.linalg.norm(points[:,:3],axis=1);keep=ci['context']&(rad>=bounds[0])&(rad<bounds[1]);xyz=np.column_stack((-points[:,1],-points[:,0],points[:,2]))
 neighborhoods=[]
 for key,title,i,total in chosen:
  center=d['xyz_Mpc'][i];q=np.flatnonzero(keep&(np.sum((xyz-center)**2,axis=1)<=(1.5*R)**2));delta=xyz[q]-center
  assert np.count_nonzero(q==d['parent_node_id'][i])==1
  # Independently recover full-ball geometry/counts at R and 1.5R.
  for ir,radius in [(1,R),(2,1.5*R)]:
   inner=(np.linalg.norm(delta,axis=1)<=radius)&(q!=d['parent_node_id'][i])
   assert inner.sum()==d[f'local_{ir}_nfull'][i]
   np.testing.assert_allclose(shape(delta[inner])[0],d[f'local_{ir}_full'][i],atol=1e-10)
  neighborhoods.append((key,title,i,total,q,delta))
 union=np.unique(np.concatenate([a[4] for a in neighborhoods]));print('selected',len(union),'neighbour rows',flush=True)
 summary=CAND/'summaries/ph006';sm=json.loads((summary/'OOF_SUMMARY_COMPLETE.json').read_text())
 ds=CAND/'dataset/P12A_DATASET_READY.json';ready=json.loads(ds.read_text())
 assert digest(summary/'OOF_SUMMARY_COMPLETE.json')==ready['phase_audit']['ph006']['summary_marker_sha256']
 assert sm['phase']=='ph006' and sm['pass'] and not sm['sealed_phase_opened']
 for name in ['parent_node_id','base_prediction','response']:
  assert digest(sm['arrays'][name])==sm['array_sha256'][name]
 parents=np.load(sm['arrays']['parent_node_id'],mmap_mode='r');lookup=np.full(len(points),-1,dtype=np.int64);lookup[parents]=np.arange(len(parents));sr=lookup[union];available=sr>=0
 cache=ROOT/'p12a_random_support_parent_cache_v2/ph006/P12A_RANDOM_SUPPORT_CACHE_READY.json';cr=json.loads(cache.read_text())
 assert cr['pass'] and not cr['sealed_phase_opened'] and cr['parent_sha256']==sm['array_sha256']['parent_node_id']
 assert digest(cr['distance_path'])==cr['distance_sha256'] and digest(cr['support_path'])==cr['support_sha256']
 support=np.load(cr['support_path'],mmap_mode='r');distance=np.load(cr['distance_path'],mmap_mode='r')
 available[available] &= support[sr[available]]
 base=np.load(sm['arrays']['base_prediction'],mmap_mode='r');response=np.load(sm['arrays']['response'])
 uu=union[available];rr=sr[available]
 ctx=np.column_stack((base[rr],response['redshift'][rr],np.log(np.maximum(response['ntilde_mpc3'][rr],1e-12)),points[uu,3],np.log1p(distance[rr]))).astype('f4')
 # Recover exact saved contexts wherever these neighbours also occur in the sample.
 sample=np.load(CAND/'dataset/ph006_selection_sample.npz');sp=sample['parent_node_id'];sl=np.full(len(points),-1,dtype=np.int64);sl[sp]=np.arange(len(sp));found=sl[uu]>=0
 np.testing.assert_allclose(ctx[found],sample['context'][sl[uu[found]]],atol=1e-6,rtol=0)
 sys.path.insert(0,'/global/u2/d/dkololgi/TNG/Illustris')
 from workflows.sbi.p12a_blind_inference import reconstruct_fmpe
 from workflows.sbi.p12_train_base_response_fmpe import sample_posterior
 cp=CAND/'posterior/fmpe_estimator.pt';assert digest(cp)==report['sha256'][str(cp)]
 posterior,ck=reconstruct_fmpe(cp,'cpu');assert ck['dataset_marker_sha256']==digest(ds)
 samples=sample_posterior(posterior,((ctx-ck['context_mean'])/ck['context_std']).astype('f4'),512,32,'cpu')
 ev=eigen(samples*np.asarray(ck['theta_std'])+np.asarray(ck['theta_mean']));classes=(ev>.2).sum(2)
 prob=np.full((len(union),4),np.nan);prob[available]=np.stack([(classes==j).mean(1) for j in range(4)],axis=1)
 # Preserve every available original evaluation probability, including all centres.
 el=np.full(len(points),-1,dtype=np.int64);el[d['parent_node_id']]=np.arange(len(t));saved=el[union]>=0
 prob[saved]=p[el[union[saved]]]
 np.testing.assert_allclose(prob[np.isfinite(prob).all(1)].sum(1),1)
 args.out.mkdir(parents=True);rows=[]
 for key,title,i,total,q,delta in neighborhoods:
  ix=np.searchsorted(union,q);e=np.asarray(native[q]);valid=ci['valid_target'][q]&np.isfinite(e).all(1);tc=np.full(len(q),-1,dtype=int);tc[valid]=(e[valid]>.2).sum(1)
  pp=prob[ix];central=q==d['parent_node_id'][i];np.testing.assert_allclose(e[central][0],d['truth_eigenvalues'][i],rtol=0,atol=1e-6);np.testing.assert_array_equal(pp[central][0],p[i])
  path=args.out/(key+'.npz');np.savez_compressed(path,parent_node_id=q,targetid=ci['targetid'][q],delta_Mpc=delta,truth=e,truth_class=tc,probability=pp,central=central,prediction_original=saved[ix])
  rows.append({'key':key,'title':title,'selection_pool_rows':total,'selection_rule':'Median parent ID in the specified geometry/truth cohort; no visual selection', 'parent_node_id':int(d['parent_node_id'][i]),'targetid':int(ci['targetid'][d['parent_node_id'][i]]),'geometry_row':i,'truth_eigenvalues':d['truth_eigenvalues'][i].tolist(),'probability':p[i].tolist(),'LPS_R':g[i].tolist(),'LPS_1p5R':d['local_2_full'][i].tolist(),'neighbors_R':int(d['local_1_nfull'][i]),'neighbors_1p5R':int(d['local_2_nfull'][i]),'projected_anisotropy':float(d['local_projected'][i]),'display_rows':len(q),'truth_available':int(valid.sum()),'prediction_available':int(np.isfinite(pp).all(1).sum()),'array_sha256':digest(path)})
 manifest={'job':os.environ['SLURM_JOB_ID'],'examples':rows,'sampling_seed':171018,'samples':512,'new_draw_rows':int(available.sum()),'reused_original_probability_rows':int(saved.sum()),'context_replay_rows':int(found.sum()),'source_sha256':{str(a):digest(a) for a in [Path(__file__),evidence/'geometry.npz',cp,ds,summary/'OOF_SUMMARY_COMPLETE.json',cache]},'scope':'Already-exposed ph006. Native host tidal labels displayed at observed galaxy positions, not a reconstructed matter field. Frozen FMPE new draws on saved encoder summaries; original evaluation probabilities retained where present. No refit. Unknown truth/predictions shown grey.'}
 (args.out/'EXAMPLES.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n');print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':main()
