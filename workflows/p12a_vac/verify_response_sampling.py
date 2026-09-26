"""Verify sparse Loa response evaluation against canonical ph006 voxel products."""
import importlib.util,json,argparse
from pathlib import Path
import numpy as np
s=importlib.util.spec_from_file_location('trial',Path(__file__).with_name('loa_trial.py'));t=importlib.util.module_from_spec(s);s.loader.exec_module(t)
from workflows.sbi.p12_prepare_base_response_dataset import sample_random_support_distance

def main(output):
 p=t.B/'training_contract_r1_random/adapters/ph006/field/adapter_manifest.json';m=t.read(p);r=t.read(m['p3_manifest']);angular=np.load(r['angular_map']);points=np.load(r['points'],mmap_mode='r');selection=t.read(r['selection_manifest']);sr=t.read(t.C/'summaries/ph006/OOF_SUMMARY_COMPLETE.json');parent=np.load(sr['arrays']['parent_node_id']);pid=parent[np.linspace(0,len(parent)-1,512,dtype=int)];ref,rs=sample_random_support_distance(m,points,pid);got=np.empty_like(ref);support=np.empty_like(rs)
 for cap,name in [(0,'SGC'),(1,'NGC')]:
  use=points[pid,3]==cap;mask=angular['support'].astype(bool)&((angular['domain']//2)==cap);angle=t.angular_boundary_distance(mask);got[use],support[use]=t.response_at(points[pid[use],:3],m['caps'][name],mask,angle,selection)
 out=dict(rows=len(pid),distance_exact=bool(np.array_equal(got,ref)),support_exact=bool(np.array_equal(support,rs)),maximum_distance_difference=float(abs(got-ref).max()),source_sha256=t.digest(Path(__file__).with_name('loa_trial.py')));t.save(output,out);print(out,flush=True)
 if not out['distance_exact'] or not out['support_exact']:raise RuntimeError('response sampling parity failed')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
