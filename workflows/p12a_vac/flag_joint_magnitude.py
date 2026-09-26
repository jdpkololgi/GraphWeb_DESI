"""All cut intersections and matched-sky z/apparent-r distributions; no repairs."""
import argparse,json,os
from pathlib import Path
import numpy as np
import fitsio,healpy as hp,h5py,hdf5plugin
from mock_selection_test import cap,digest,save
REPO=Path(__file__).resolve().parents[2]
PREV=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_mock_test_ph006_20260926_v1')
BASE=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2')
REAL=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
ZE=np.linspace(.15,.55,41);ME=np.linspace(12,20.2,83)
FLAGS=['BGS_BRIGHT_bit','GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED','COADD_FIBERSTATUS_zero','DELTACHI2_gt40_keep_GALAXY','uniform_r_lt19p5','nominal_photsys_magnitude']

def intersections(pattern):
    codes=np.arange(256)
    return np.stack([pattern[..., (codes & c)==c].sum(-1) for c in codes],axis=-1)

def histogram(z,m,c,k):
    return np.array([np.histogram2d(z[k&(c==i)],m[k&(c==i)],bins=[ZE,ME])[0] for i in (0,1)],dtype='i8')

def run(root):
    if (root/'COMPLETE.json').exists():raise FileExistsError('completed run exists')
    common=np.load(PREV/'common.npy');patterns=np.zeros((2,40,256),dtype='i8');real_hist=np.zeros((2,40,82),dtype='i8');strict_hist=real_hist.copy();fibre={};fail={};total=0
    cols=['TARGETID','RA','DEC','Z_not4clus','ZWARN','SPECTYPE','DELTACHI2','BGS_TARGET','GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED','COADD_FIBERSTATUS','FLUX_R','MW_TRANSMISSION_R','PHOTSYS']
    source_before=REAL.stat()
    with fitsio.FITS(REAL) as f:
      for start in range(0,f[1].get_nrows(),250000):
        d=f[1].read(rows=np.arange(start,min(start+250000,f[1].get_nrows())),columns=cols)
        z=d['Z_not4clus'];base=(z>=.15)&(z<.55)&(d['ZWARN']==0)&(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
        pix=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True);d=d[base&common[pix]];z=d['Z_not4clus'];c=cap(d['RA'],d['DEC']);total+=len(d)
        valid=(d['FLUX_R']>0)&(d['MW_TRANSMISSION_R']>0);m=np.full(len(d),np.nan);m[valid]=22.5-2.5*np.log10(d['FLUX_R'][valid]/d['MW_TRANSMISSION_R'][valid]);north=np.char.strip(d['PHOTSYS'].astype('U'))=='N'
        conditions=[(d['BGS_TARGET']&2)!=0,d['GOODHARDLOC'].astype(bool),d['GOODPRI'].astype(bool),d['LOCATION_ASSIGNED'].astype(bool),d['COADD_FIBERSTATUS']==0,d['DELTACHI2']>40,valid&(m<19.5),valid&(m>=12)&(m<np.where(north,19.54,19.5))]
        code=np.zeros(len(d),dtype='i8')
        for j,k in enumerate(conditions):code|=k.astype('i8')<<j
        zi=np.searchsorted(ZE,z,side='right')-1
        patterns+=np.bincount((c*40+zi)*256+code,minlength=2*40*256).reshape(2,40,256)
        real_hist+=histogram(z,m,c,np.ones(len(d),bool));strict_hist+=histogram(z,m,c,code==255)
        for value,n in zip(*np.unique(d['COADD_FIBERSTATUS'],return_counts=True)):fibre[str(int(value))]=fibre.get(str(int(value)),0)+int(n)
        for j,k in enumerate(conditions):fail[FLAGS[j]]=fail.get(FLAGS[j],0)+int((~k).sum())
        if start%2000000==0:print('Loa',start,flush=True)
    assert (source_before.st_size,source_before.st_mtime_ns)==(REAL.stat().st_size,REAL.stat().st_mtime_ns)
    combo=intersections(patterns);assert int(combo[...,0].sum())==total
    # True parent magnitudes/redshifts for successful AND unassigned mock targets.
    parent=BASE/'forFA6_nomask.fits';parent_hist=np.zeros_like(real_hist)
    with fitsio.FITS(parent) as f:
      n=f[1].get_nrows();zlookup=np.full(n+1,np.nan);mlookup=np.full(n+1,np.nan)
      for start in range(0,n,500000):
        d=f[1].read(rows=np.arange(start,min(start+500000,n)),columns=['TARGETID','BGS_TARGET','RA','DEC','RSDZ','R_MAG_APP'])
        ids=d['TARGETID'];assert np.array_equal(ids,np.arange(start+1,start+len(d)+1))
        zlookup[ids]=d['RSDZ'];mlookup[ids]=d['R_MAG_APP'];c=cap(d['RA'],d['DEC']);k=((d['BGS_TARGET']&2)!=0)&common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
        parent_hist+=histogram(d['RSDZ'],d['R_MAG_APP'],c,k)
    mock=BASE/'altmtl6/loa-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.h5';full_hist=np.zeros_like(real_hist);obs_hist=np.zeros_like(real_hist);mock_flags={};matched=0
    with h5py.File(mock) as f:
      t=f['LSS'];n=len(t['TARGETID'])
      for start in range(0,n,250000):
        sl=slice(start,min(start+250000,n));d={k:t[k][sl] for k in ['TARGETID','RA','DEC','ZWARN','BGS_TARGET','GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED']};ids=d['TARGETID'];assert ((ids>0)&(ids<len(zlookup))).all();z=zlookup[ids];m=mlookup[ids];assert np.isfinite(z).all() and np.isfinite(m).all();matched+=len(ids)
        c=cap(d['RA'],d['DEC']);k=common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)];good=d['ZWARN']==0
        full_hist+=histogram(z,m,c,k);obs_hist+=histogram(z,m,c,k&good)
        for name in ['GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED']:
          mock_flags[name]=mock_flags.get(name,0)+int((k&good&(z>=.15)&(z<.55)&~d[name].astype(bool)).sum())
    expected=json.loads((REPO/'docs/evidence/p12a_closure_20260926/PRODUCT_CROSSWALK.json').read_text())['products']
    assert np.array_equal(combo[...,0],expected['real_loa']['counts'])
    assert np.array_equal(obs_hist.sum(-1),expected['mock_loa']['counts'])
    assert np.all(obs_hist<=full_hist) and np.all(full_hist<=parent_hist)
    xyz=hp.pix2ang(256,np.flatnonzero(common),lonlat=True);cc=cap(*xyz);area=np.bincount(cc,minlength=2)*hp.nside2pixarea(256,degrees=True)
    rows=[]
    for j in range(256):
      sh=combo[...,j].reshape(2,4,10).sum(-1);den=obs_hist.sum(-1).reshape(2,4,10).sum(-1)
      rows.append(dict(code=j,cuts=[FLAGS[k] for k in range(8) if j&(1<<k)],counts=sh.tolist(),loa_to_mock=(sh/den).tolist(),removed=int(total-combo[...,j].sum())))
    np.savez_compressed(root/'HISTOGRAMS.npz',z_edges=ZE,r_edges=ME,patterns=patterns,combinations=combo,loa=real_hist,loa_strict=strict_hist,mock_parent=parent_hist,mock_full=full_hist,mock_observed=obs_hist,area_deg2=area)
    decoded={}
    try:
      from desispec.maskbits import fibermask
      decoded={k:fibermask.names(int(k)) for k in fibre}
    except ImportError:decoded={'status':'desispec unavailable; integer flags preserved'}
    save(root/'COMPLETE.json',dict(flags=FLAGS,combinations=rows,baseline_rows=total,single_cut_failures=fail,fibre_status_counts=fibre,fibre_status_names=decoded,mock_common_success_flag_failures=mock_flags,mock_joined_rows=matched,area_deg2=area.tolist(),
       sources={'loa':str(REAL),'mock_full':str(mock),'parent':str(parent)},common_mask_sha256=digest(PREV/'common.npy'),script_sha256=digest(__file__),histogram_sha256=digest(root/'HISTOGRAMS.npz'),job=os.environ['SLURM_JOB_ID'],
       notes=['All 256 conjunctions use the existing ZWARN0/GALAXY/DELTACHI2>=25 baseline; GALAXY is never silently dropped.', 'Nominal magnitude conjunction is a diagnostic and can reject legitimate SGA recovery rows.', 'Mock lacks observed DELTACHI2, SPECTYPE and COADD_FIBERSTATUS counterparts; no synthetic values assigned.', 'Loa retention is conditional on already successful baseline objects, not absolute completeness.', 'Mock retention recovers true parent z for unassigned targets; includes angular vetoes when parent is denominator.', 'Apparent-magnitude passbands and estimators remain unproven equivalent.'],validation={'previous_counts_exact':True,'mock_retention_nested':True}))
    print('COMPLETE',total,fail,'strict ratios',rows[-1]['loa_to_mock'],flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True);run(a.root)
