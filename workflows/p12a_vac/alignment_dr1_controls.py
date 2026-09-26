"""DR1/Loa random-supported footprint and same-TARGETID transition controls."""
import argparse,os
from pathlib import Path
import fitsio
import healpy as hp
import numpy as np
from mock_selection_test import cap,digest,save
from flag_joint_magnitude import REAL
from alignment_parent_screen import MASK,ZE
IRON=Path('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/BGS_BRIGHT_full_HPmapcut.dat.fits')
COLS=['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE','FLUX_R','MW_TRANSMISSION_R','LOCATION_ASSIGNED','PHOTSYS']

def read_data(path):
    parts=[]
    with fitsio.FITS(path) as f:
        for start in range(0,f[1].get_nrows(),250000):
            d=f[1].read(rows=np.arange(start,min(start+250000,f[1].get_nrows())),columns=COLS)
            o=np.empty(len(d),dtype=[('id','i8'),('z','f8'),('r','f8'),('zw','?'),('gal','?'),('dc','?'),('assigned','?'),('cap','i1'),('north','?'),('pix','i8')])
            o['id']=d['TARGETID'];o['z']=d['Z_not4clus'];o['r']=np.nan
            good=(d['FLUX_R']>0)&(d['MW_TRANSMISSION_R']>0)
            o['r'][good]=22.5-2.5*np.log10(d['FLUX_R'][good]/d['MW_TRANSMISSION_R'][good])
            o['zw']=d['ZWARN']==0;o['gal']=np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY'
            o['dc']=d['DELTACHI2']>=25;o['assigned']=d['LOCATION_ASSIGNED'].astype(bool)
            o['north']=np.char.strip(d['PHOTSYS'].astype('U'))=='N'
            o['cap']=cap(d['RA'],d['DEC']);o['pix']=hp.ang2pix(512,d['RA'],d['DEC'],lonlat=True,nest=True)
            parts.append(o)
    o=np.concatenate(parts);o.sort(order='id');assert np.all(np.diff(o['id'])>0)
    print('indexed',path.parent,len(o),flush=True);return o

def quality(t):return t['zw']&t['gal']&t['dc']&(t['z']>=.15)&(t['z']<.55)
def count(t,k):return np.array([np.histogram(t['z'][k&(t['cap']==i)],ZE)[0] for i in [0,1]])
def random_map(path):
    h=np.zeros(hp.nside2npix(512),dtype='i8')
    with fitsio.FITS(path) as f:
        n=f[1].get_nrows()
        for start in range(0,n,500000):
            d=f[1].read(rows=np.arange(start,min(start+500000,n)),columns=['RA','DEC'])
            pix=hp.ang2pix(512,d['RA'],d['DEC'],lonlat=True,nest=True)
            h+=np.bincount(pix,minlength=len(h))
    assert h.sum()==n;print('randoms',path.parent,n,flush=True);return h

def run(root):
    root.mkdir(parents=True,exist_ok=True)
    if (root/'COMPLETE.json').exists():raise FileExistsError(root)
    paths={'iron':IRON,'loa':REAL};rans={k:p.with_name('BGS_BRIGHT_0_full_HPmapcut.ran.fits') for k,p in paths.items()}
    before={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in [*paths.values(),*rans.values()]}
    maps={k:random_map(p) for k,p in rans.items()}
    old=hp.reorder(np.load(MASK),r2n=True).astype(bool)
    masks={'original':np.repeat(old,4),'random256':np.repeat(old & (maps['iron'].reshape(-1,4).sum(-1)>0)&(maps['loa'].reshape(-1,4).sum(-1)>0),4),
           'random512':np.repeat(old,4)&(maps['iron']>0)&(maps['loa']>0),
           'random512_interior':np.repeat(old,4)&(maps['iron']>=16)&(maps['loa']>=16)}
    a=read_data(IRON);b=read_data(REAL);qa=quality(a);qb=quality(b)
    hist={};foot={}
    for name,m in masks.items():
        foot[name]={'occupied_pixel_area_deg2':float(m.sum()*hp.nside2pixarea(512,degrees=True)),
                    'random_rows':{k:int(v[m].sum()) for k,v in maps.items()}}
        for key,t,q in [('iron',a,qa),('loa',b,qb)]:
            sky=m[t['pix']];hist[key+'_'+name]=count(t,q&sky)
            hist[key+'_'+name+'_r19p5']=count(t,q&sky&(t['r']>=12)&(t['r']<19.5))
    ai,bi=np.intersect1d(a['id'],b['id'],assume_unique=True,return_indices=True)[1:]
    aa=a[ai];bb=b[bi];qaa=qa[ai];qbb=qb[bi]
    matched_a=np.zeros(len(a),bool);matched_a[ai]=True
    matched_b=np.zeros(len(b),bool);matched_b[bi]=True
    m=masks['random512'];sky=m[aa['pix']]&m[bb['pix']]
    hist['matched_iron']=count(aa,qaa&sky);hist['matched_loa']=count(bb,qbb&sky)
    hist['iron_only_ids']=count(a,qa&m[a['pix']]&~matched_a)
    hist['loa_only_ids']=count(b,qb&m[b['pix']]&~matched_b)
    # Exact histogram decomposition allows migration of jointly accepted redshifts.
    hist['gained_at_matched_ids']=count(bb,~qaa&qbb&sky)
    hist['lost_at_matched_ids']=count(aa,qaa&~qbb&sky)
    hist['both_iron_z']=count(aa,qaa&qbb&sky);hist['both_loa_z']=count(bb,qaa&qbb&sky)
    reasons={};remainder=~qaa&qbb&sky
    for name,cond in [('not_assigned',~aa['assigned']),('ZWARN_nonzero',~aa['zw']),('SPECTYPE_not_GALAXY',~aa['gal']),('DELTACHI2_below25',~aa['dc']),('redshift_outside',~((aa['z']>=.15)&(aa['z']<.55)))]:
        take=remainder&cond;hist['gain_reason_'+name]=count(bb,take);reasons[name]=int(take.sum());remainder &= ~take
    assert not remainder.any()
    # Categorical success accounting, plus exact conservation on matched sky.
    assert np.array_equal(hist['matched_loa']-hist['matched_iron'],hist['gained_at_matched_ids']-hist['lost_at_matched_ids']+hist['both_loa_z']-hist['both_iron_z'])
    both=qaa&qbb&sky
    stats={'matched_ids':len(ai),'both_science':int(both.sum()),'both_science_abs_dr_gt_001':int((both&(abs(aa['r']-bb['r'])>.01)).sum()),
           'both_science_relative_dz_gt_0005':int((both&(abs(aa['z']-bb['z'])/(1+aa['z'])>.005)).sum()),
           'matched_pixel_changes':int((aa['pix']!=bb['pix']).sum()),'gain_reasons_exclusive':reasons}
    assert before=={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in [*paths.values(),*rans.values()]}
    np.savez_compressed(root/'HISTOGRAMS.npz',**hist,z_edges=ZE)
    save(root/'COMPLETE.json',dict(counts={k:v.reshape(2,4,10).sum(-1).tolist() for k,v in hist.items()},footprints=foot,stats=stats,
        sources=before,job=os.environ.get('SLURM_JOB_ID'),script_sha256=digest(__file__),histogram_sha256=digest(root/'HISTOGRAMS.npz'),mask_sha256=digest(MASK),
        helpers={n:digest(Path(__file__).with_name(n)) for n in ['mock_selection_test.py','flag_joint_magnitude.py','alignment_parent_screen.py']},
        contract='Full BGS_BRIGHT iron v1.5 vs Loa v2.1; Z_not4clus, ZWARN0, DCHI>=25, GALAXY, .15<=z<.55. No weights. Random index0, HEALPix NEST512 with256 aggregation; interior >=16 random rows per pixel in each release. Not exact angular masks or area estimates. Gains classified at matched IDs using Loa z; missing IDs are catalogue membership differences, not automatically new observations.',
        limitations=['One random realization; interior threshold is a sensitivity control, not survey selection.','No claim of independently certified release/blinding state.','Conditional reobservation diagnostics do not measure redshifts of all failures.']))
    print('COMPLETE',stats,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
