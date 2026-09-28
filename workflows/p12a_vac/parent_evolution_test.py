"""Fixed-parent evolution intervention, diagnostic only; exposed ph000."""
import os,json,hashlib,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
import healpy as hp
from alignment_parent_screen import BASE,MASK
R=Path(__file__).resolve().parents[2]
O=R/'docs/evidence/p12a_evolution_test_20260928';O.mkdir(exist_ok=True)
sky=np.load(MASK);ze=np.linspace(.15,.55,5)
counts={};metadata={};stride=20
for version in ['v0.1','v1']:
 path=BASE/version/'z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
 before=[path.stat().st_size,path.stat().st_mtime_ns]
 variants=['baseline'] if version=='v0.1' else ['baseline','global_Q078','colour_Q','Q0']
 h={k:np.zeros((2,4),dtype=int) for k in variants};near=np.zeros((2,4),dtype=int)
 cols=['RA','DEC','Z','R_MAG_APP','G_R_REST']
 with fits.open(path,memmap=True) as f:
  for start in range(0,len(f[1].data),2000000):
   d={c:np.array(f[1].data[c][start:start+2000000:stride]) for c in cols}
   z=d['Z'];r=d['R_MAG_APP'];sel=(z>=.15)&(z<.55)&sky[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
   ps=d['DEC']>32.375
   for v in variants:
    q={'baseline':.67,'global_Q078':.78,'Q0':0}.get(v,np.where(d['G_R_REST']>=.75,.23,1.59))
    rr=r if v=='baseline' else r-(q-.67)*(z-.1)
    for i in [0,1]:h[v][i]+=np.histogram(z[sel&(ps==i)&(rr>=12)&(rr<19.5)],ze)[0]
   for i in [0,1]:near[i]+=np.histogram(z[sel&(ps==i)&(r>=19.5)&(r<19.92)],ze)[0]
   print(version,start,flush=True)
 assert before==[path.stat().st_size,path.stat().st_mtime_ns]
 counts[version]={k:v.tolist() for k,v in h.items()}
 metadata[version]={'path':str(path),'stat':before,'available_fainter_rows_19p5_19p92':near.tolist()}
ratios={k:(np.array(v)/np.array(counts['v1']['baseline'])).tolist() for k,v in counts['v1'].items()}
result={'counts':counts,'v1_over_baseline':ratios,'sources':metadata,'z_edges':ze.tolist(),'hemispheres':['photometric S','photometric N'],'stride':stride,'job':os.environ.get('SLURM_JOB_ID'),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitations':'Fixed existing parent, colour and K; no LF/HOD refit. Colour split uses stored G_R_REST and assumes published reference convention. No guarantee missing faint parent rows can be recovered. No observation pipeline; no correction applied. Deterministic stride20, no significance estimate.'}
(O/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(ratios))
