"""Descriptive colour comparison at common magnitude limit and matched fine z."""
import json,hashlib
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_targeting_alignment_20260925.json'
def main():
 d=json.loads(P.read_text());lo=np.array(d['catalogues']['loa']['fine_gr_hist_uniform_r19p5']);mo=np.array(d['catalogues']['ph006']['fine_gr_hist_uniform_r19p5']);edges=np.array(d['gr_edges']);cen=(edges[1:]+edges[:-1])/2;out=[]
 for cap in [0,1]:
  for shell in range(4):
   l=lo[cap,shell*10:(shell+1)*10];m=mo[cap,shell*10:(shell+1)*10];assert np.all(m.sum(1)>0);mw=m*(l.sum(1)/m.sum(1))[:,None];lh=l.sum(0);mh=mw.sum(0);med=lambda x:float(np.interp(.5*x.sum(),np.cumsum(x),cen))
   out.append(dict(cap=cap,shell=shell,desi_median_gr=med(lh),mock_z_reweighted_median_gr=med(mh),median_difference=med(lh)-med(mh),desi_rows=int(lh.sum()),mock_rows=int(m.sum())))
 d['colour_uniform_r19p5_matched_fine_z']=out;d['colour_summary_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();P.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
