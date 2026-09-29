"""Exact-ID join of supplied Loa CIGALE catalogue to the unchanged posterior VAC."""
import argparse,json,os
from pathlib import Path
import numpy as np
import fitsio
from property_environment import VAC,digest
SOURCE=Path('/global/cfs/cdirs/desicollab/users/zouhu/vac/dr2/dr2_galaxy_sedfitting_v1.0.fits')
PROPS=[f'{x}_CG_{fit}' for fit in ['15','5'] for x in ['MASS','MASSERR','SFR','SFRERR','AGE','AGEERR','AV','AVERR','FAGN','FAGNERR']]
COLS=['TARGETID','TARGET_RA','TARGET_DEC','Z','SPECTYPE','SURVEY','PROGRAM','TNAME','SNR_MED','FLUX_SCALE']+PROPS

def main(out):
    out.mkdir(parents=True,exist_ok=False);source_stat=SOURCE.stat();before=digest(VAC)
    d=fitsio.read(VAC,columns=['TARGETID','RA','DEC','Z','SUPPORTED','QUALITY','P_WEB'])
    ids=d['TARGETID'];order=np.argsort(ids);sid=ids[order];assert np.all(np.diff(sid)>0)
    pieces=[];total=0
    with fitsio.FITS(SOURCE) as f:
        h=f[1];n=h.get_nrows();(out/'SOURCE_HEADER.txt').write_text(str(h.read_header()))
        for lo in range(0,n,500000):
            c=h.read(rows=np.arange(lo,min(lo+500000,n)),columns=COLS)
            ix=np.searchsorted(sid,c['TARGETID']);ok=(ix<len(sid))&(sid[np.minimum(ix,len(sid)-1)]==c['TARGETID'])
            if ok.any():pieces.append((c[ok],order[ix[ok]],lo+np.flatnonzero(ok)));total+=int(ok.sum())
            if lo%5000000==0:print('Scanned',lo,'matched rows',total,flush=True)
    c=np.concatenate([v[0] for v in pieces]);vr=np.concatenate([v[1] for v in pieces]);sr=np.concatenate([v[2] for v in pieces]);del pieces
    dz=abs(c['Z']-d['Z'][vr])/(1+d['Z'][vr]);dra=(c['TARGET_RA']-d['RA'][vr]+180)%360-180
    sep=np.hypot(dra*np.cos(np.deg2rad(d['DEC'][vr])),c['TARGET_DEC']-d['DEC'][vr])*3600
    consistent=np.isfinite(dz)&np.isfinite(sep)&(dz<.001)&(sep<1)
    unique_before=len(np.unique(vr));bad_ids=c['TARGETID'][~consistent]
    c=c[consistent];vr=vr[consistent];sr=sr[consistent];dz=dz[consistent];sep=sep[consistent]
    idx=np.argsort(vr,kind='stable');c=c[idx];vr=vr[idx];sr=sr[idx];dz=dz[idx];sep=sep[idx]
    unique,first,count=np.unique(vr,return_index=True,return_counts=True)
    # Ambiguous repeated observations are excluded rather than choosing an SED arbitrarily.
    keep=count==1;ambiguous=unique[~keep];ii=first[keep];c=c[ii];vr=vr[ii];sr=sr[ii];dz=dz[ii];sep=sep[ii]
    dtype=list(d.dtype.descr)+[('CIGALE_MATCHED','?'),('CIGALE_AMBIGUOUS','?'),('CIGALE_SOURCE_ROW','i8'),('CIGALE_Z','f8'),('CIGALE_SEP_ARCSEC','f8'),('CIGALE_DZ_NORM','f8')]+[(f'CIGALE_{x}',c.dtype[x]) for x in ['SPECTYPE','SURVEY','PROGRAM','TNAME','SNR_MED','FLUX_SCALE']]+[(x,'f8') for x in PROPS]
    joined=np.zeros(len(d),dtype=dtype)
    for name in d.dtype.names:joined[name]=d[name]
    for name in joined.dtype.names:
        if name not in d.dtype.names and joined.dtype[name].kind=='f':joined[name]=np.nan
    joined['CIGALE_SOURCE_ROW']=-1;joined['CIGALE_MATCHED'][vr]=True;joined['CIGALE_AMBIGUOUS'][ambiguous]=True;joined['CIGALE_SOURCE_ROW'][vr]=sr
    joined['CIGALE_Z'][vr]=c['Z'];joined['CIGALE_SEP_ARCSEC'][vr]=sep;joined['CIGALE_DZ_NORM'][vr]=dz
    for name in ['SPECTYPE','SURVEY','PROGRAM','TNAME','SNR_MED','FLUX_SCALE']:joined[f'CIGALE_{name}'][vr]=c[name]
    for name in PROPS:joined[name][vr]=c[name]
    dest=out/'DESI_LOA_P12A_CIGALE_JOIN.fits';fitsio.write(dest,joined,header={'PROVIS':True,'CGSOURCE':str(SOURCE),'CGMETHOD':'Unique TARGETID, sep<1 arcsec, abs(dz)/(1+z)<0.001'},clobber=False)
    # Check output identity and exact frozen posterior retention, not only row counts.
    check=fitsio.read(dest,columns=['TARGETID','P_WEB','SUPPORTED','QUALITY'])
    for name in check.dtype.names:assert np.array_equal(check[name],d[name],equal_nan=True)
    report={'source':str(SOURCE),'source_rows':n,'source_sha256':digest(SOURCE),'source_bytes':source_stat.st_size,
            'vac':str(VAC),'vac_sha256':before,'vac_rows':len(d),'id_matched_source_rows':total,'unique_vac_ids_matched_before_coordinate_check':unique_before,
            'coordinate_redshift_rejected_rows':len(bad_ids),'rejected_targetids':bad_ids.tolist(),'ambiguous_duplicate_ids':int(len(ambiguous)),
            'accepted_unique_matches':len(c),'unmatched_or_excluded_vac_rows':len(d)-len(c),'max_sep_arcsec':float(sep.max()),'max_dz_normalized':float(dz.max()),
            'columns':COLS,'joined_path':str(dest),'joined_sha256':digest(dest),'source_stat_unchanged':SOURCE.stat().st_mtime_ns==source_stat.st_mtime_ns,
            'vac_unchanged':digest(VAC)==before,'posterior_columns_exactly_preserved':True,'job_id':os.environ.get('SLURM_JOB_ID'),'script_sha256':digest(__file__)}
    assert report['vac_unchanged'] and report['source_stat_unchanged']
    (out/'JOIN.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rejected_targetids'},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)
