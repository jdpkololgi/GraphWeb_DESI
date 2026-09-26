"""Recover mock target redshifts by TARGETID, including unassigned sentinel rows."""
import numpy as np,fitsio
from pathlib import Path
import loa_trial as t
from selection_audit import ROOT,MOCK
receipt=t.read(t.B/'ph006/catalogues/observed/ph006_bgs_bright_full_observed_with_tweb.fits.complete.json')
parent=receipt['annotated_parent']['path']
with fitsio.FITS(parent) as f:
 z=np.empty(f[1].get_nrows(),dtype='f8')
 for start in range(0,len(z),250000):
  d=f[1][start:min(start+250000,len(z))]
  assert np.array_equal(d['TARGETID'],np.arange(start+1,start+len(d)+1))
  z[start:start+len(d)]=d['Z']
counts={k:np.zeros(8,dtype='i8') for k in ['full_targets','zwarn0','assigned','goodhardware','goodpriority','valid_source_z','valid_truez']}
with fitsio.FITS(MOCK) as f:
 for start in range(0,f[1].get_nrows(),250000):
  d=f[1][start:min(start+250000,f[1].get_nrows())];ids=d['TARGETID'];assert np.all((ids>0)&(ids<=len(z)));tz=z[ids-1];sh=np.searchsorted([.15,.25,.35,.45,.55],tz,side='right')-1;ok=(sh>=0)&(sh<4);cap=t.galactic_cap(d['RA'],d['DEC']);b=cap*4+np.clip(sh,0,3)
  for name,mask in dict(full_targets=np.ones(len(d),bool),zwarn0=d['ZWARN']==0,assigned=d['LOCATION_ASSIGNED'].astype(bool),goodhardware=d['GOODHARDLOC'].astype(bool),goodpriority=d['GOODPRI'].astype(bool),valid_source_z=(d['Z_not4clus']>=.15)&(d['Z_not4clus']<.55),valid_truez=(d['TRUEZ']>=.15)&(d['TRUEZ']<.55)).items():counts[name]+=np.bincount(b[ok&mask],minlength=8)
report=dict(parent=parent,mock=MOCK,counts={k:v.tolist() for k,v in counts.items()},note='All full-table targets binned by parent observed/RSD redshift through TARGETID; unlike full-table TRUEZ and Z_not4clus this includes unassigned sentinel rows. No label columns used.')
t.save(ROOT/'MOCK_ASSIGNMENT_AUDIT.json',report);print(report)
