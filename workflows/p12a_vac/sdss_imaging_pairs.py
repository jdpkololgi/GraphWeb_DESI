"""Bounded direct SDSS imaging crossmatch of the paired-SED pilot targets."""
from pathlib import Path
import os,json,hashlib
import numpy as np,fitsio
from scipy.spatial import cKDTree
from same_galaxy_passbands import OUT,CAT,quant
BASE=Path('/global/cfs/cdirs/cosmo/data/sdss/dr16/eboss/sweeps/dr13_final')
def xyz(ra,dec):
 r=np.deg2rad(ra);d=np.deg2rad(dec);return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])
def main():
 s=np.load(OUT/'SAME_GALAXY_SAMPLE.npz');rows=s['catalogue_rows']
 with fitsio.FITS(CAT) as f:m=f['METADATA'].read(rows=rows,columns=['TARGETID','RA','DEC','Z','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R'])
 assert np.array_equal(m['TARGETID'],s['targetid'])
 ix=fitsio.read(BASE/'datasweep-index-gal.fits',ext=1);tree=cKDTree(xyz(ix['RA'],ix['DEC']));near=tree.query_ball_point(xyz(m['RA'],m['DEC']),2*np.sin(np.deg2rad(.3)/2));files={}
 for i,ns in enumerate(near):
  for j in ns:
   row=ix[j]
   if row['NPRIMARY']<=0:continue
   key=(str(row['RERUN']).strip(),int(row['RUN']),int(row['CAMCOL']));v=files.setdefault(key,{'targets':set(),'fields':set()});v['targets'].add(i);v['fields'].add(j)
 chosen=sorted(files,key=lambda k:(-len(files[k]['targets']),k))[:16];records=[];receipts=[]
 for rerun,run,cam in chosen:
  key=(rerun,run,cam);p=BASE/rerun/f'calibObj-{run:06d}-{cam}-gal.fits.gz';ranges=[np.arange(int(ix[j]['ISTART']),int(ix[j]['IEND'])+1) for j in files[key]['fields']];rr=np.unique(np.concatenate(ranges))
  with fitsio.FITS(p) as f:
   assert rr.max()<f[1].get_nrows();d=f[1].read(rows=rr,columns=['RA','DEC','RESOLVE_STATUS','THING_ID','EXTINCTION','MODELFLUX','MODELFLUX_IVAR','PETROFLUX','PETROFLUX_IVAR','FLAGS'])
  d=d[(d['RESOLVE_STATUS']&256)!=0]
  if not len(d):continue
  q=cKDTree(xyz(d['RA'],d['DEC']));dist,at=q.query(xyz(m['RA'],m['DEC']),k=2);arc=np.rad2deg(2*np.arcsin(np.minimum(dist/2,1)))*3600
  for i in np.flatnonzero((arc[:,0]<1)&(arc[:,1]>=1)):
   x=d[at[i,0]];pf=x['PETROFLUX'][1:3];mf=x['MODELFLUX'][1:3]
   if np.any(pf<=0) or np.any(mf<=0):continue
   legacy=22.5-2.5*np.log10([m[i]['FLUX_G'],m[i]['FLUX_R']]);pm=22.5-2.5*np.log10(pf)-x['EXTINCTION'][1:3];mm=22.5-2.5*np.log10(mf)-x['EXTINCTION'][1:3]
   sn=pf*np.sqrt(np.maximum(x['PETROFLUX_IVAR'][1:3],0))
   records.append(dict(targetid=int(m[i]['TARGETID']),z=float(m[i]['Z']),sep_arcsec=float(arc[i,0]),thing_id=int(x['THING_ID']),path=str(p),petro_minus_legacy_r=float(pm[1]-legacy[1]),model_minus_legacy_r=float(mm[1]-legacy[1]),petro_minus_model_r=float(pm[1]-mm[1]),petro_minus_legacy_gr=float(pm[0]-pm[1]-legacy[0]+legacy[1]),petro_snr_gr=sn.tolist(),flags_gr=[int(v) for v in x['FLAGS'][1:3]]))
  receipts.append(dict(path=str(p),rows_read=len(rr),targets_near_fields=len(files[key]['targets'])));print(run,cam,len(records),flush=True)
 # Reject duplicate matches across loaded files rather than silently choose.
 from collections import Counter
 counts=Counter(r['targetid'] for r in records);unique=[r for r in records if counts[r['targetid']]==1];summary=[]
 for lo,hi in zip([.15,.25,.35,.45],[.25,.35,.45,.55]):
  a=[r for r in unique if lo<=r['z']<hi and min(r['petro_snr_gr'])>=5]
  summary.append(dict(z=[lo,hi],n=len(a),stats={k:quant(np.array([r[k] for r in a])) for k in ['petro_minus_legacy_r','model_minus_legacy_r','petro_minus_model_r','petro_minus_legacy_gr'] if a}))
 out=dict(job=os.environ.get('SLURM_JOB_ID'),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),candidate_files=len(files),loaded_files=receipts,total_pilot_targets=len(m),raw_matches=len(records),unique_matches=len(unique),duplicate_targets=sum(v>1 for v in counts.values()),summary=summary,records=unique,method='At most 16 sweep files with most pilot targets within .3deg of field center. SURVEY_PRIMARY bit256; unique match within1arcsec and no second within1arcsec; positive g/r flux; summaries require Petrosian SNR>=5 both bands. SDSS stored extinction versus already de-reddened FastSpecFit Legacy flux. No full photometric clean-bit selection or aperture correction; pilot footprint incomplete, unmatched fraction not completeness.')
 (OUT/'SDSS_IMAGING_PAIRS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
