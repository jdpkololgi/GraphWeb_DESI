"""Isolated, exposed-ph006 selection experiment; never a production/quality model."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import fitsio
import healpy as hp
import numpy as np

REPO = Path(__file__).resolve().parents[2]
BASE = Path('/pscratch/sd/d/dkololgi/abacus/p10_multiphase')
CANARY = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
RAW = Path('/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v0.1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph006.fits')
EDGES = np.linspace(.15, .55, 41)

def digest(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''): h.update(b)
    return h.hexdigest()

def save(p, a):
    p = Path(p)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(a, indent=2, allow_nan=False) + '\n')
    tmp.replace(p)

def cap(ra, dec):
    # Same rotation as the P3/P12 cap convention, with no dependence on colour.
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    return (SkyCoord(ra=ra*u.deg, dec=dec*u.deg).galactic.b.deg > 0).astype('i1')

def counts(z, c, keep):
    return np.array([np.histogram(z[keep & (c == k)], EDGES)[0] for k in (0, 1)])

def quality(d):
    good = d['ZWARN'] == 0
    return {'quality25': good & (d['DELTACHI2'] >= 25) & (np.char.strip(d['SPECTYPE'].astype('U')) == 'GALAXY'),
            'quality40': good & (d['DELTACHI2'] > 40)}

def trial_mask(z, m, delta, q):
    return m + delta - q * (z - .1) < 19.5

def source(p):
    s = Path(p).stat()
    return dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns)

def prepare(root):
    if (root/'PREPARED.json').exists(): raise FileExistsError('prepared root already exists')
    common = np.load(CANARY/'galaxy_angular_support.npy')
    mask_sources = [source(CANARY/'galaxy_angular_support.npy')]
    for i in range(2, 7):
        p = BASE/f'ph{i:03d}/p3_fields/angular_support_nside256.npz'
        common &= np.load(p)['support']; mask_sources.append(source(p))
    domain = np.load(CANARY/'random_angular.npz')['domain']
    np.save(root/'common.npy', common)
    raw_cols = ['RA','DEC','Z','R_MAG_APP','R_MAG_ABS','G_R_OBS','HALO_MASS','CEN','FILE_NUM','HALO_INDEX','BOX_INDEX','IN_Y1','IN_Y5']
    parts = []
    before = source(RAW)
    with fitsio.FITS(RAW) as f:
        n = f[1].get_nrows()
        for start in range(0, n, 500000):
            d = f[1][start:min(start+500000,n)]
            pix = hp.ang2pix(256, d['RA'], d['DEC'], lonlat=True)
            ok = common[pix] & (d['Z'] >= .15) & (d['Z'] < .55)
            at = np.flatnonzero(ok)
            a = np.empty(len(at), dtype=[('RAW_ROW','i8'),('CAP','i1'),('PHOTSYS','i1')] + [(k,d.dtype[k]) for k in raw_cols])
            a['RAW_ROW'] = start + at
            a['CAP'] = cap(d['RA'][ok],d['DEC'][ok])
            a['PHOTSYS'] = domain[pix[ok]] % 2
            for k in raw_cols: a[k] = d[k][ok]
            parts.append(a)
            if start % 10000000 == 0: print('raw',start,n,flush=True)
    raw = np.concatenate(parts); del parts
    assert source(RAW) == before
    assert np.all(np.diff(raw['RAW_ROW']) > 0)
    np.save(root/'raw_common.npy', raw)
    loa = Path(json.loads((CANARY/'INPUTS_READY.json').read_text())['source']['path'])
    before_loa = source(loa); real_parts=[]; stages={}
    cols = ['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE','PHOTSYS','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R']
    with fitsio.FITS(loa) as f:
        for start in range(0,f[1].get_nrows(),250000):
            d = f[1].read(rows=np.arange(start,min(start+250000,f[1].get_nrows())),columns=cols)
            z=d['Z_not4clus']; pix=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)
            ok=common[pix] & np.isfinite(z) & (z>=.15) & (z<.55)
            d=d[ok]; z=z[ok]; c=cap(d['RA'],d['DEC']); masks=quality(d)
            valid=(d['FLUX_G']>0)&(d['FLUX_R']>0)&(d['MW_TRANSMISSION_G']>0)&(d['MW_TRANSMISSION_R']>0)
            r=np.full(len(d),np.nan);colour=r.copy()
            r[valid]=22.5-2.5*np.log10(d['FLUX_R'][valid]/d['MW_TRANSMISSION_R'][valid])
            colour[valid]=-2.5*np.log10(d['FLUX_G'][valid]/d['MW_TRANSMISSION_G'][valid])-(r[valid]-22.5)
            for name,k in masks.items():
                for suffix,cut in [('all_targeted',np.ones(len(d),bool)),('r19p5',valid&(r>=12)&(r<19.5))]:
                    key=name+'_'+suffix;stages.setdefault(key,np.zeros((2,40),dtype='i8'));stages[key]+=counts(z,c,k&cut)
            k=masks['quality25']&valid&(r>=12)&(r<19.5)
            a=np.empty(k.sum(),dtype=[('TARGETID','i8'),('Z','f8'),('R','f8'),('COLOUR','f8'),('CAP','i1'),('PHOTSYS','S1')])
            for key,val in [('TARGETID',d['TARGETID']),('Z',z),('R',r),('COLOUR',colour),('CAP',c),('PHOTSYS',d['PHOTSYS'])]:a[key]=val[k]
            real_parts.append(a)
    assert source(loa)==before_loa
    real=np.concatenate(real_parts); assert len(np.unique(real['TARGETID']))==len(real)
    np.save(root/'loa_common.npy',real)
    raw_counts=counts(raw['Z'],raw['CAP'],raw['R_MAG_APP']<19.5)
    previous=json.loads((REPO/'docs/figures/p12a_selection_20260925/UPSTREAM.json').read_text())
    assert np.array_equal(raw_counts,np.array(previous['counts']['raw_rlt19.5'])[:,:40])
    old=json.loads((REPO/'docs/evidence/p12a_closure_20260926/PRODUCT_CROSSWALK.json').read_text())
    assert np.array_equal(stages['quality25_all_targeted'],old['products']['real_loa']['counts'])
    save(root/'PREPARED.json',dict(raw=before,loa=before_loa,masks=mask_sources,raw_rows=len(raw),loa_rows=len(real),
        stages={k:v.tolist() for k,v in stages.items()},raw_counts=raw_counts.tolist(),baseline_counts_replay=True,
        script_sha256=digest(__file__),job=os.environ.get('SLURM_JOB_ID'),arrays={p.name:digest(p) for p in [root/'raw_common.npy',root/'loa_common.npy',root/'common.npy']}))
    print('PREPARED',len(raw),len(real),flush=True)

def shell_counts(a): return a.reshape(2,4,10).sum(-1)

def sweep(root, config):
    prep=json.loads((root/'PREPARED.json').read_text()); raw=np.load(root/'raw_common.npy',mmap_mode='r');real=np.load(root/'loa_common.npy',mmap_mode='r')
    for name,sha in prep['arrays'].items():
        if digest(root/name)!=sha:raise ValueError('prepared input hash changed: '+name)
    for f in ['SWEEP.json','trial_parent.fits']:
        if (root/f).exists():raise FileExistsError(f)
    target=np.array(prep['stages']['quality25_r19p5']);base=np.array(prep['raw_counts'])
    obs=np.array(json.loads((REPO/'docs/evidence/p12a_closure_20260926/MOCK_RETENTION.json').read_text())['stage']['loa_zwarn0'])
    retention=np.divide(obs,base); assert np.all((retention>=0)&(retention<=1))
    target_shell=shell_counts(target);rows=[]
    def evaluate(delta,q,stage):
        k=trial_mask(raw['Z'],raw['R_MAG_APP'],delta,q)
        h=counts(raw['Z'],raw['CAP'],k);pred=h*retention;ratio=shell_counts(pred)/target_shell
        score=float(np.mean(np.log(ratio[0])**2))
        row=dict(delta_m=float(delta),extra_q=float(q),stage=stage,sgc_log_mse=score,ngc_log_mse=float(np.mean(np.log(ratio[1])**2)),
            raw_counts=h.tolist(),fixed_retention_prediction=pred.tolist(),predicted_to_loa_shell_ratio=ratio.tolist(),
            raw_parent_deficient_shells=(shell_counts(h)<target_shell).tolist())
        rows.append(row);return row
    for delta in config['delta_m_grid']:
        for q in config['extra_q_grid']:evaluate(delta,q,'coarse')
    best=min(rows,key=lambda r:r['sgc_log_mse']);d0=best['delta_m'];q0=best['extra_q']
    for delta in np.linspace(max(-.10,d0-.05),min(.10,d0+.05),5):
        for q in np.linspace(max(-.2,q0-.2),min(.8,q0+.2),5):
            if any(abs(r['delta_m']-delta)<1e-8 and abs(r['extra_q']-q)<1e-8 for r in rows):continue
            evaluate(delta,q,'refine')
    best=min(rows,key=lambda r:r['sgc_log_mse']);summaries={}
    for name,delta,q in [('baseline',0.,0.),('screening_best',best['delta_m'],best['extra_q'])]:
        m=raw['R_MAG_APP']+delta-q*(raw['Z']-.1);k=trial_mask(raw['Z'],raw['R_MAG_APP'],delta,q);ss=[]
        for c in (0,1):
            for j in range(4):
                lo=.15+.1*j;hi=lo+.1;mk=k&(raw['CAP']==c)&(raw['Z']>=lo)&(raw['Z']<hi)
                rk=(real['CAP']==c)&(real['Z']>=lo)&(real['Z']<hi)
                ss.append(dict(cap=c,z=[lo,hi],mock_count=int(mk.sum()),real_count=int(rk.sum()),
                    mock_colour_median=float(np.median(raw['G_R_OBS'][mk])),real_colour_median=float(np.median(real['COLOUR'][rk])),
                    mock_r_median=float(np.median(m[mk])),real_r_median=float(np.median(real['R'][rk]))))
        summaries[name]=ss
    best_keep=trial_mask(raw['Z'],raw['R_MAG_APP'],best['delta_m'],best['extra_q'])
    data=raw[best_keep];extra=np.empty(len(data),dtype=data.dtype.descr+[('R_MAG_TRIAL','f8'),('BGS_TARGET_TRIAL','i8')])
    for name in data.dtype.names:extra[name]=data[name]
    extra['R_MAG_TRIAL']=data['R_MAG_APP']+best['delta_m']-best['extra_q']*(data['Z']-.1);extra['BGS_TARGET_TRIAL']=2
    fitsio.write(root/'trial_parent.fits',extra,header={'PHASE':'ph006','TESTONLY':True,'D_M':best['delta_m'],'Q_EXTRA':best['extra_q'],'FARERUN':False})
    result=dict(config=config,config_sha256=digest(REPO/'configs/p12a_mock_test/ph006_20260926.json'),script_sha256=digest(__file__),
        trials=rows,best=best,population=summaries,retention=retention.tolist(),target_counts=target.tolist(),
        trial_parent=dict(path=str(root/'trial_parent.fits'),rows=len(extra),sha256=digest(root/'trial_parent.fits')),
        limitations=['Fixed baseline shell retention is a screening approximation, not rerun assignment or a spectral quality model.',
                      'Mock and real passbands/flux estimators are not established equivalent.',
                      'SGC tuning and NGC check are exposed-data diagnostics, not independent validation.',
                      'Stored absolute-magnitude evolution convention unresolved; extra_q is not physical LF Q.'],
        science_ready=False,job=os.environ.get('SLURM_JOB_ID'))
    save(root/'SWEEP.json',result); print('BEST',json.dumps(best),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','sweep']);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    config=json.loads((REPO/'configs/p12a_mock_test/ph006_20260926.json').read_text())
    if config['phase']!='ph006':raise PermissionError('only exposed ph006 is allowed')
    a.root.mkdir(parents=True,exist_ok=True)
    if a.mode=='prepare':prepare(a.root)
    else:sweep(a.root,config)
if __name__=='__main__':main()
