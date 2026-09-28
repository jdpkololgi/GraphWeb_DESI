"""Validate prepared BGS target types, magnitude split and finite subpriorities."""
import json,sys,hashlib
from pathlib import Path
import fitsio,numpy as np,healpy as hp
from mock_selection_test import digest,cap
from alignment_parent_screen import MASK
sky=np.load(MASK);hist=np.zeros((2,40),dtype=np.int64);edges=np.linspace(.15,.55,41)
p=Path(sys.argv[1]);counts={'bright':0,'faint':0,'faint_hip':0}
with fitsio.FITS(p) as f:
 n=f[1].get_nrows();last=0
 for start in range(0,n,250000):
  d=f[1].read(rows=np.arange(start,min(start+250000,n)),columns=['TARGETID','BGS_TARGET','R_MAG_APP','DEC','SUBPRIORITY','ZWARN','RA','RSDZ'])
  assert d['TARGETID'][0]==last+1 and np.all(np.diff(d['TARGETID'])==1);last=int(d['TARGETID'][-1])
  assert np.all(np.isfinite(d['SUBPRIORITY'])&(d['SUBPRIORITY']>=0)&(d['SUBPRIORITY']<1))
  b=d['BGS_TARGET']==2;faint=d['BGS_TARGET']==1;hip=d['BGS_TARGET']==9
  assert np.all(b|faint|hip)
  cut=np.where(d['DEC']>32.375,19.54,19.5)
  assert np.all(d['R_MAG_APP'][b]<cut[b])
  assert np.all((d['R_MAG_APP'][~b]>=cut[~b])&(d['R_MAG_APP'][~b]<=20.175))
  assert np.all(d['ZWARN']==0)
  common=b&(d['R_MAG_APP']>=12)&(d['R_MAG_APP']<19.5)&(d['RSDZ']>=.15)&(d['RSDZ']<.55)
  common &= sky[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
  c=cap(d['RA'],d['DEC'])
  for i in [0,1]:hist[i]+=np.histogram(d['RSDZ'][common&(c==i)],edges)[0]
  for key,mask in [('bright',b),('faint',faint),('faint_hip',hip)]:counts[key]+=int(mask.sum())
assert sum(counts.values())==n
(p.parent/'PREPARATION_READY.json').write_text(json.dumps({'target':str(p),'rows':n,'types':counts,'sha256':digest(p),'placeholder_zwarn':True,'qualification':False,'common_uniform_bright_counts':hist.tolist(),'common_uniform_bright_shell_counts':hist.reshape(2,4,10).sum(-1).tolist(),'definition':'BRIGHT r in [12,19.5), RSDZ in [.15,.55), common NSIDE256 sky; after tile footprint and imaging mask, before assignment/veto/redshift success'},indent=2)+'\n')
