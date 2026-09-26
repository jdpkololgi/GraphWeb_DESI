"""Internal DA3 vs Loa ID/quality crosswalk; no release substitution."""
import argparse
import os
from pathlib import Path
import fitsio
import numpy as np
import healpy as hp
from mock_selection_test import cap,digest,save
from flag_joint_magnitude import REAL
from alignment_parent_screen import MASK,ZE

DA3=Path('/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/LSScats/v0/BGS_BRIGHT_full_HPmapcut.dat.fits')
COLS=['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE','FLUX_R','MW_TRANSMISSION_R']

def compact(d,mask):
    z=d['Z_not4clus'];m=np.full(len(d),np.nan)
    good=(d['FLUX_R']>0)&(d['MW_TRANSMISSION_R']>0)
    m[good]=22.5-2.5*np.log10(d['FLUX_R'][good]/d['MW_TRANSMISSION_R'][good])
    q=(d['ZWARN']==0)&(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')&(z>=.15)&(z<.55)
    out=np.empty(len(d),dtype=[('id','i8'),('z','f8'),('r','f8'),('q','?'),('sky','?'),('cap','i1')])
    out['id']=d['TARGETID'];out['z']=z;out['r']=m;out['q']=q
    out['sky']=mask[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)];out['cap']=cap(d['RA'],d['DEC'])
    return out

def count(t,k):
    return np.array([np.histogram(t['z'][k&(t['cap']==i)],bins=ZE)[0] for i in (0,1)])

def run(out):
    if (out/'COMPLETE.json').exists():raise FileExistsError(out)
    mask=np.load(MASK);pieces=[];meta={}
    with fitsio.FITS(REAL) as f:
        n=f[1].get_nrows()
        for start in range(0,n,250000):pieces.append(compact(f[1].read(rows=np.arange(start,min(start+250000,n)),columns=COLS),mask))
    a=np.concatenate(pieces);del pieces
    a.sort(order='id');assert np.all(np.diff(a['id'])>0)
    print('Loa indexed',len(a),flush=True)
    counts={'loa_common_sky':count(a,a['q']&a['sky'])}
    for name in ['da3_common_sky','matched_loa_common_sky','matched_da3_common_sky','da3_unmatched_common_sky']:
        counts[name]=np.zeros((2,40),dtype='i8')
    stats=dict(matched_ids=0,matched_both_science=0,loa_science_lost_at_matched_ids=0,da3_science_gained_at_matched_ids=0,matched_both_science_dr_gt_001=0,matched_both_science_relative_dz_gt_0005=0)
    b_ids=[]
    with fitsio.FITS(DA3) as f:
        n=f[1].get_nrows();h=f[1].read_header()
        meta['da3_header_non_schema']={k:h[k] for k in h.keys() if not k.startswith(('TTYPE','TFORM','TUNIT','TDISP'))}
        for start in range(0,n,250000):
            b=compact(f[1].read(rows=np.arange(start,min(start+250000,n)),columns=COLS),mask);b_ids.append(b['id'].copy())
            counts['da3_common_sky']+=count(b,b['q']&b['sky'])
            ix=np.searchsorted(a['id'],b['id']);match=(ix<len(a));ix=np.minimum(ix,len(a)-1);match&=a['id'][ix]==b['id']
            counts['da3_unmatched_common_sky']+=count(b,~match&b['q']&b['sky'])
            aa=a[ix[match]];bb=b[match];sky=aa['sky']&bb['sky']
            counts['matched_loa_common_sky']+=count(aa,aa['q']&sky)
            counts['matched_da3_common_sky']+=count(bb,bb['q']&sky)
            both=aa['q']&bb['q']&sky
            stats['matched_ids']+=int(match.sum())
            stats['matched_both_science']+=int(both.sum())
            stats['loa_science_lost_at_matched_ids']+=int((aa['q']&~bb['q']&sky).sum())
            stats['da3_science_gained_at_matched_ids']+=int((~aa['q']&bb['q']&sky).sum())
            stats['matched_both_science_dr_gt_001']+=int((both&(abs(aa['r']-bb['r'])>.01)).sum())
            stats['matched_both_science_relative_dz_gt_0005']+=int((both&(abs(aa['z']-bb['z'])/(1+aa['z'])>.005)).sum())
    ids=np.concatenate(b_ids);assert len(np.unique(ids))==len(ids)
    np.savez_compressed(out/'HISTOGRAMS.npz',**counts,z_edges=ZE)
    save(out/'COMPLETE.json',dict(counts={k:v.reshape(2,4,10).sum(-1).tolist() for k,v in counts.items()},stats=stats,metadata=meta,
        sources={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in [REAL,DA3]},job=os.environ.get('SLURM_JOB_ID'),
        common_mask_sha256=digest(MASK),script_sha256=digest(__file__),histogram_sha256=digest(out/'HISTOGRAMS.npz'),
        notes=['Same full-catalogue Z_not4clus quality definition applied. No FKP weighting or clustering positions used.',
               'Blinding/version completeness not independently certified; differences are delivered-column diagnostics, not cosmological evolution.',
               'Matched-ID quality transitions are conditional on membership in both delivered full catalogues and original common occupied sky.',
               'Unmatched IDs can reflect masks/catalogue construction as well as new observations; do not label all as newly observed galaxies.',
               'No environmental model execution or change of primary Loa target.']))
    print('COMPLETE release',stats,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True);run(a.root)
