"""Observed/rest-frame Loa colour audit; requires approved compute execution."""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import socket
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / 'loa_mass_complete_20261016'
sys.path.insert(0, str(PREVIOUS))
# Import the frozen estimator under a distinct name: this script is also run.py.
import importlib.util
_spec = importlib.util.spec_from_file_location('previous_mass_analysis', PREVIOUS/'run.py')
_previous = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_previous)
controlled = _previous.controlled
SOURCE, PHOT = _previous.OUT, _previous.PHOT
BASE = SOURCE.parent
OUT = BASE/'loa_rest_colour_20261016_v1'
FSF = Path('/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs')
WINDOWS = [('z015020', .15, .20, 10.60, 534070), ('z020022', .20, .22, 10.69, 232024),
           ('z020030', .20, .30, 10.96, 611197)]
CUTS = (.70, .75, .80)
ENV = ['Void', 'Wall', 'Filament', 'Cluster']
SPECCOLS = ['TARGETID', 'ABSMAG01_SDSS_G', 'ABSMAG01_SDSS_R',
            'ABSMAG01_IVAR_SDSS_G', 'ABSMAG01_IVAR_SDSS_R', 'RCHI2_CONT']
METACOLS = ['TARGETID', 'RA', 'DEC', 'Z', 'Z_RR', 'ZWARN_RR', 'DELTACHI2']

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()

def write(name, obj):
    (OUT/name).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')

def masks(c, cut=.75):
    return {'red': np.isfinite(c)&(c>cut)&(c<1.5), 'blue': np.isfinite(c)&(c<cut)}

def category(c, cut=.75):
    v=np.full(len(c),2,dtype=np.int8)
    mm=masks(c,cut);v[mm['blue']]=0;v[mm['red']]=1
    return v

def index_match(ids, query):
    at=np.searchsorted(ids,query);safe=np.minimum(at,len(ids)-1)
    return safe, (at<len(ids))&(ids[safe]==query)

def sep_arcsec(ra1,dec1,ra2,dec2):
    dr=(ra1-ra2+180)%360-180
    return 3600*np.hypot(dr*np.cos(np.deg2rad((dec1+dec2)/2)),dec1-dec2)

def percentile(values):
    q=float(np.quantile(values,.95));v=np.sort(values);pos=.95*(len(v)-1);j=int(pos)
    assert abs(q-(v[j]+(pos-j)*(v[min(j+1,len(v)-1)]-v[j])))<1e-10
    return q

def counts(a, colour='GR_REST01', cut=.75):
    mm=masks(a[colour],cut);hard=a['P_WEB'].argmax(1)
    return dict(rows=len(a), colours={k:int(v.sum()) for k,v in mm.items()},
        hard={k:np.bincount(hard[v],minlength=4).tolist() for k,v in mm.items()},
        excluded=int((~(mm['red']|mm['blue'])).sum()))

def parent_join():
    import fitsio
    b0=np.load(SOURCE/'bright_parent.npz')
    keep=(b0['z']>=.15)&(b0['z']<.30)
    order=np.argsort(b0['TARGETID'][keep]);b={k:b0[k][keep][order] for k in b0.files}
    ids=b['TARGETID'];n=len(ids);assert (np.diff(ids)>0).all()
    obs=np.full(n,np.nan);ra=np.full(n,np.nan);dec=np.full(n,np.nan);seen=np.zeros(n,np.int16)
    cols=['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R']
    with fitsio.FITS(PHOT) as f:
        for start in range(0,f[1].get_nrows(),300000):
            x=f[1].read(rows=np.arange(start,min(start+300000,f[1].get_nrows())),columns=cols)
            at,ok=index_match(ids,x['TARGETID']);ok&=(x['ZWARN']==0)&(x['DELTACHI2']>40)&(x['Z_not4clus']>=.15)&(x['Z_not4clus']<.30)
            x=x[ok];at=at[ok];np.add.at(seen,at,1)
            assert np.allclose(x['Z_not4clus'],b['z'][at],rtol=0,atol=1e-12)
            ra[at]=x['RA'];dec[at]=x['DEC']
            valid=(x['FLUX_G']>0)&(x['FLUX_R']>0)&(x['MW_TRANSMISSION_G']>0)&(x['MW_TRANSMISSION_R']>0)
            with np.errstate(divide='ignore',invalid='ignore'):
                col=-2.5*np.log10((x['FLUX_G']/x['MW_TRANSMISSION_G'])/(x['FLUX_R']/x['MW_TRANSMISSION_R']))
            obs[at]=np.where(valid,col,np.nan)
    assert (seen==1).all()
    dtype=[('TARGETID','i8'),('FSF_FILE','i2'),('FSF_ROW','i8'),('FSF_NMATCH','i2'),
           ('FSF_Z','f8'),('FSF_SEP_ARCSEC','f8'),('FSF_DZ_NORM','f8'),
           ('REST_VALID','?'),('GR_REST01','f8'),('GR_OBS','f8')]
    dtype += [(c,'f8') for c in SPECCOLS if c!='TARGETID']
    matched=np.zeros(n,dtype=dtype);matched['TARGETID']=ids;matched['FSF_FILE']=-1;matched['FSF_ROW']=-1
    for c in matched.dtype.names:
        if matched[c].dtype.kind=='f': matched[c]=np.nan
    matched['GR_OBS']=obs
    files=sorted(FSF.glob('fastspec-loa-main-bright-nside1-hp*.fits'));assert len(files)==12
    provenance=[]
    for fi,path in enumerate(files):
        h=hashlib.sha256();nmatch=0
        with fitsio.FITS(path) as f:
            assert f['METADATA'].get_nrows()==f['SPECPHOT'].get_nrows()
            for start in range(0,f['METADATA'].get_nrows(),300000):
                meta=f['METADATA'].read(rows=np.arange(start,min(start+300000,f['METADATA'].get_nrows())),columns=METACOLS)
                at,ok=index_match(ids,meta['TARGETID']); rows=start+np.flatnonzero(ok)
                if not len(rows):continue
                spec=f['SPECPHOT'].read(rows=rows,columns=SPECCOLS);meta=meta[ok];at=at[ok]
                assert np.array_equal(meta['TARGETID'],spec['TARGETID'])
                h.update(meta.tobytes());h.update(spec.tobytes());np.add.at(matched['FSF_NMATCH'],at,1)
                matched['FSF_FILE'][at]=fi;matched['FSF_ROW'][at]=rows;matched['FSF_Z'][at]=meta['Z']
                sep=sep_arcsec(ra[at],dec[at],meta['RA'],meta['DEC']);dz=np.abs(meta['Z']-b['z'][at])/(1+b['z'][at])
                matched['FSF_SEP_ARCSEC'][at]=sep;matched['FSF_DZ_NORM'][at]=dz
                for c in SPECCOLS:
                    if c!='TARGETID':matched[c][at]=spec[c]
                gr=spec['ABSMAG01_SDSS_G'].astype(float)-spec['ABSMAG01_SDSS_R'].astype(float)
                valid=np.isfinite(gr)&np.isfinite(spec['ABSMAG01_IVAR_SDSS_G'])&np.isfinite(spec['ABSMAG01_IVAR_SDSS_R'])
                valid&=(spec['ABSMAG01_IVAR_SDSS_G']>0)&(spec['ABSMAG01_IVAR_SDSS_R']>0)
                valid&=(spec['ABSMAG01_SDSS_G']!=0)&(spec['ABSMAG01_SDSS_R']!=0)&np.isfinite(spec['RCHI2_CONT'])&(spec['RCHI2_CONT']>0)
                valid&=(meta['ZWARN_RR']==0)&(meta['DELTACHI2']>40)&(sep<=1)&(dz<=.001)
                valid&=(np.abs(meta['Z_RR']-b['z'][at])/(1+b['z'][at])<=.001)
                matched['REST_VALID'][at]=valid;matched['GR_REST01'][at]=gr;nmatch+=len(at)
        provenance.append(dict(index=fi,path=str(path),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,matched_rows=nmatch,matched_input_columns_sha256=h.hexdigest()))
        print('FSF',path.name,nmatch,flush=True)
    matched['REST_VALID']&=matched['FSF_NMATCH']==1
    matched['GR_REST01'][~matched['REST_VALID']]=np.nan
    np.save(OUT/'parent_rest_join.npy',matched,allow_pickle=False)
    # Fresh independent source reads for a stratified sample across all shards.
    validated=0
    for fi,path in enumerate(files):
        ii=np.flatnonzero((matched['FSF_FILE']==fi)&matched['REST_VALID'])
        if not len(ii):continue
        ii=ii[np.linspace(0,len(ii)-1,min(32,len(ii)),dtype=int)]
        fresh=fitsio.read(path,ext='SPECPHOT',rows=matched['FSF_ROW'][ii],columns=SPECCOLS)
        assert np.array_equal(fresh['TARGETID'],matched['TARGETID'][ii])
        for c in SPECCOLS:
            if c!='TARGETID':assert np.array_equal(fresh[c],matched[c][ii]),c
        assert np.array_equal(fresh['ABSMAG01_SDSS_G'].astype(float)-fresh['ABSMAG01_SDSS_R'].astype(float),matched['GR_REST01'][ii])
        validated+=len(ii)
    write('join_audit.json',dict(parent_rows=n,valid=int(matched['REST_VALID'].sum()),unmatched=int((matched['FSF_NMATCH']==0).sum()),
        duplicate=int((matched['FSF_NMATCH']>1).sum()),invalid_single=int(((matched['FSF_NMATCH']==1)&~matched['REST_VALID']).sum()),
        source_row_checks=validated,files=provenance,source_parent_sha256=sha(SOURCE/'bright_parent.npz'),
        definition='ABSMAG01_SDSS_G - ABSMAG01_SDSS_R; SDSS bands shifted to zref=0.1, not zref=0; no internal dust correction',
        quality='Unique main-bright match, both magnitude IVAR>0, finite nonzero magnitudes, RCHI2_CONT>0 finite, ZWARN_RR=0, DELTACHI2>40, sep<=1arcsec, normalized z and z_RR agreement<=0.001'))
    return b,matched

def completeness(b,matched):
    old=json.loads((SOURCE/'mass_limits.json').read_text());audit=[];offset=b['r']-b['rlim']
    for prev in old:
        if prev['zlo']<.15 or prev['zhi']>.30:continue
        sel=(b['z']>=prev['zlo'])&(b['z']<prev['zhi']);limits=[prev['logmass_limit']];pops={}
        for cut in CUTS:
            for group,cm in masks(matched['GR_REST01'],cut).items():
                ii=np.flatnonzero(sel&cm);assert len(ii)>=100,(prev['zlo'],cut,group,len(ii))
                faint=ii[offset[ii]>=np.quantile(offset[ii],.8)]
                lim=percentile(b['mass'][faint]+.4*offset[faint]);limits.append(lim)
                pops[f'{cut:.2f}_{group}']=dict(rows=len(ii),faint_rows=len(faint),p95=lim)
        audit.append(dict(zlo=prev['zlo'],zhi=prev['zhi'],base_limit=prev['logmass_limit'],rest_limit=max(limits),
            parent_rows=int(sel.sum()),rest_valid=int((sel&matched['REST_VALID']).sum()),populations=pops))
    assert len(audit)==15
    write('rest_completeness.json',audit)
    return audit

def transfer(a, good):
    rec=dict(total=len(a),valid=int(good.sum()),fraction=float(good.mean()))
    for name,labels in [('hard',a['P_WEB'].argmax(1)),('zbin',np.floor((a['Z']-.15+1e-9)/.01).astype(int)),
                        ('massbin',np.floor((np.log10(a['MASS_CG_15'])-6)/.1).astype(int))]:
        rec[name]=[{ 'bin':int(j),'total':int((labels==j).sum()),'valid':int(((labels==j)&good).sum())} for j in np.unique(labels)]
    return rec

def prepare():
    import fitsio
    from numpy.lib.recfunctions import append_fields
    OUT.mkdir(exist_ok=False)
    b,matched=parent_join();curve=completeness(b,matched)
    chunks=[]
    with fitsio.FITS(SOURCE/'LOA_P12A_MASS_COMPLETE_zdependent.fits') as f:
        for start in range(0,f[1].get_nrows(),200000):
            a=f[1].read(rows=np.arange(start,min(start+200000,f[1].get_nrows())))
            chunks.append(a[(a['Z']>=.15)&(a['Z']<.30)])
    a=np.concatenate(chunks);del chunks
    at,ok=index_match(matched['TARGETID'],a['TARGETID']);assert ok.all()
    assert np.allclose(a['GR_OBS'],matched['GR_OBS'][at],rtol=0,atol=1e-12,equal_nan=True)
    add=[c for c in matched.dtype.names if c not in ['TARGETID','GR_OBS']]
    a=append_fields(a,add,[matched[c][at] for c in add],usemask=False)
    windows={};mass=np.log10(a['MASS_CG_15'])
    for key,lo,hi,oldlim,expected in WINDOWS:
        bins=[v for v in curve if v['zlo']>=lo and v['zhi']<=hi]
        restlim=float(np.ceil(max(v['rest_limit'] for v in bins)*100)/100)
        assert restlim>=max(v['base_limit'] for v in bins)
        pairlim=max(restlim,oldlim)
        zw=(a['Z']>=lo)&(a['Z']<hi)
        old=a[zw&(mass>=oldlim)];assert len(old)==expected
        fitsio.write(OUT/f'observed_original_{key}.fits',old)
        primary=a[zw&(mass>=restlim)&a['REST_VALID']]
        paired=a[zw&(mass>=pairlim)&a['REST_VALID']]
        for label,data in [('rest_primary',primary),('paired',paired)]:
            fitsio.write(OUT/f'{label}_{key}.fits',data,header={'MCFIX':restlim if label=='rest_primary' else pairlim,'MCZLO':lo,'MCZHI':hi,'RESTZ':.1,'COLCUT':.75})
            ix=np.linspace(0,len(data)-1,min(257,len(data)),dtype=int)
            fresh=fitsio.read(OUT/f'{label}_{key}.fits',rows=ix)
            for col in data.dtype.names:
                assert np.array_equal(fresh[col],data[col][ix],equal_nan=True) if data[col].dtype.kind in 'fc' else np.array_equal(fresh[col],data[col][ix])
        transition=np.bincount(category(paired['GR_OBS'])*3+category(paired['GR_REST01']),minlength=9).reshape(3,3)
        oldvalid=old[old['REST_VALID']]
        oldtrans=np.bincount(category(oldvalid['GR_OBS'])*3+category(oldvalid['GR_REST01']),minlength=9).reshape(3,3)
        sf=(paired['SFR_CG_15']>=0)&np.isfinite(paired['SFR_CG_15'])
        fq={}
        for cut in CUTS:
            fq[str(cut)]={}
            for name,cm in masks(primary['GR_REST01'],cut).items():
                use=cm&np.isfinite(primary['SFR_CG_15'])&(primary['SFR_CG_15']>=0)
                fq[str(cut)][name]=dict(rows=int(cm.sum()),valid_sfr=int(use.sum()),quenched_fraction=float((primary['SFR_CG_15'][use]/primary['MASS_CG_15'][use]<1e-11).mean()))
        w=dict(zlo=lo,zhi=hi,old_limit=oldlim,rest_limit=restlim,paired_limit=pairlim,
               observed_original=counts(old,'GR_OBS'),rest_primary=counts(primary),paired_observed=counts(paired,'GR_OBS'),paired_rest=counts(paired),
               match_coverage=transfer(old,old['REST_VALID']),original_mass_same_rows_observed=counts(oldvalid,'GR_OBS'),
               original_mass_same_rows_rest=counts(oldvalid),transition_old_mass=oldtrans.tolist(),transition_paired=transition.tolist(),
               transition_order=['blue','red','outside'],colour_vs_quenching=fq,
               median_observed_minus_rest=float(np.median(oldvalid['GR_OBS']-oldvalid['GR_REST01'])))
        n=w['rest_primary'];w['rest_count_gate']=min(n['colours'].values())>=20000 and min(min(v) for v in n['hard'].values())>=100
        windows[key]=w;print('WINDOW',key,json.dumps(w),flush=True)
    write('windows.json',windows)
    # Hold the mass floor fixed across all redshift bins to distinguish the
    # colour-frame effect from the changing window-specific mass selections.
    common_limit=max(w['paired_limit'] for w in windows.values())
    diagnostic=[]
    for bb in curve:
        use=(a['Z']>=bb['zlo'])&(a['Z']<bb['zhi'])&(mass>=common_limit)&a['REST_VALID']
        d=a[use];rec=dict(zlo=bb['zlo'],zhi=bb['zhi'],rows=len(d))
        for frame,col in [('observed','GR_OBS'),('rest','GR_REST01')]:
            mm=masks(d[col]);rec[frame]={k:int(v.sum()) for k,v in mm.items()}
            rec[frame]['blue_fraction_all']=float(mm['blue'].mean()) if len(d) else None
        diagnostic.append(rec)
    write('fixed_mass_colour_fractions.json',dict(logmass=common_limit,bins=diagnostic,
        denominator='All valid matched galaxies above the same mass floor in each bin, including those outside red/blue cuts'))
    write('prepare_manifest.json',dict(job=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),code_sha256=sha(Path(__file__)),
        input_vac_sha256=sha(SOURCE/'LOA_P12A_MASS_COMPLETE_zdependent.fits'),outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file()}))

# Original estimator's bootstrap and common-stratum support retained exactly.
def raw(y,w,region):
    if not len(y) or (w.sum(0)<=0).any():return None
    _,reg=np.unique(region,return_inverse=True);nr=reg.max()+1
    c=np.array([np.bincount(reg,weights=w[:,j],minlength=nr) for j in range(4)])
    s=np.array([np.bincount(reg,weights=w[:,j]*y,minlength=nr) for j in range(4)])
    mean=s.sum(1)/c.sum(1);boot=[];rng=np.random.default_rng(20261016)
    for _ in range(256):
        mult=np.bincount(rng.integers(nr,size=nr),minlength=nr);den=c@mult
        if (den>0).all():boot.append((s@mult)/den)
    if len(boot)<205:return None
    boot=np.array(boot);assert np.allclose(mean,w.T@y/w.sum(0))
    return dict(mean=mean.tolist(),ci95=np.quantile(boot,[.025,.975],axis=0).tolist(),contrast=float(mean[-1]-mean[0]),
                contrast_ci95=np.quantile(boot[:,-1]-boot[:,0],[.025,.975]).tolist(),bootstrap_valid=len(boot))

def analyse(data,colour,cut=.75,properties=None):
    import healpy as hp
    results={}
    for group,mask in masks(data[colour],cut).items():
        a=data[mask];m=np.log10(a['MASS_CG_15']);z=a['Z'];sf=a['SFR_CG_15'];valid=np.isfinite(sf)&(sf>=0)
        p=a['P_WEB'].astype(float);assert np.isfinite(p).all() and (p>=0).all() and np.allclose(p.sum(1),1,atol=1e-5)
        hard=np.eye(4)[p.argmax(1)];region=hp.ang2pix(8,a['RA'],a['DEC'],lonlat=True)
        with np.errstate(divide='ignore',invalid='ignore'):
            ss=np.where(valid,np.log10(sf/a['MASS_CG_15']),np.nan)
            props={'quenched':np.where(valid,(ss<-11).astype(float),np.nan),'logmass':m,'logssfr':np.maximum(ss,-14),
                   'logsfr':np.where(valid,np.maximum(np.log10(sf),-14),np.nan),'colour':a[colour],
                   'logage':np.log10(np.where(a['AGE_CG_15']>0,a['AGE_CG_15']/1000,np.nan)),
                   'av':np.where(a['AV_CG_15']>=0,a['AV_CG_15'],np.nan),
                   'fagn':np.where((a['FAGN_CG_15']>=0)&(a['FAGN_CG_15']<=1),a['FAGN_CG_15'],np.nan)}
        block={}
        for prop,y in props.items():
            if properties is not None and prop not in properties:continue
            ii=np.flatnonzero(np.isfinite(y));rec={'rows':len(ii),'raw':{method:raw(y[ii],w[ii],region[ii]) for method,w in [('soft',p),('hard',hard)]},'controlled':{}}
            for width in [.1,.05]:
                zb=np.floor((z[ii]-.2+1e-9)/.01).astype(int);cells=zb if prop=='logmass' else np.floor((m[ii]-6)/width).astype(int)*100+zb
                for mode,weights in [('soft',{'soft':p[ii]}),('joint',{'soft':p[ii],'hard':hard[ii]})]:
                    r=controlled(y[ii],weights,cells,region[ii]) if len(ii) else {'status':'unavailable: empty population'}
                    if r['status']=='measured':
                        kept=np.isin(cells,r['common_cell_ids']);r['matched_z_range']=[float(z[ii][kept].min()),float(z[ii][kept].max())]
                        for method,w in weights.items():
                            check=np.zeros(4)
                            for cell,mix in zip(r['common_cell_ids'],r['reference_weights']):
                                use=cells==cell;check+=mix*(w[use].T@y[ii][use])/w[use].sum(0)
                            assert np.allclose(check,r['methods'][method]['mean'],atol=1e-11)
                    rec['controlled'][str(width)+'_'+mode]=r
            block[prop]=rec
        results[group]=block
    return results

def worker(task):
    import fitsio
    key,variant,cut=task
    paired=variant.startswith('paired');file=OUT/f"{'paired' if paired else 'rest_primary'}_{key}.fits"
    cols=['MASS_CG_15','SFR_CG_15','Z','P_WEB','RA','DEC','GR_OBS','GR_REST01','AGE_CG_15','AV_CG_15','FAGN_CG_15']
    a=fitsio.read(file,columns=cols);colour='GR_OBS' if variant=='paired_observed' else 'GR_REST01'
    result=analyse(a,colour,cut,properties=None if cut==.75 else ['quenched','logssfr','colour'])
    name=f'relationships_{variant}_{key}_cut{round(cut*100):02d}.json';write(name,result)
    print('ANALYSED',name,flush=True);return name

def run_analysis():
    tasks=[(key,variant,.75) for key,*_ in WINDOWS for variant in ['rest_primary','paired_observed','paired_rest']]
    tasks += [(key,'rest_primary',cut) for key,*_ in WINDOWS for cut in [.7,.8]]
    with ProcessPoolExecutor(max_workers=6) as pool:files=list(pool.map(worker,tasks))
    write('analysis_manifest.json',dict(job=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),code_sha256=sha(Path(__file__)),
        estimator_sha256=sha(PREVIOUS/'run.py'),outputs={f:sha(OUT/f) for f in files},
        validation='All raw and controlled mean vectors independently reconstructed; 256 attempted sky bootstraps; unsupported controls explicitly unavailable',
        limitations='Conditional on selection, colour cuts and inferred environment; excludes SED, colour-boundary, matching/selection and simulation-transfer systematics; no causal or calibration claim'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['prepare','analyse','all'],default='all');args=ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):raise RuntimeError('Approved compute allocation required')
    if args.stage in ['prepare','all']:prepare()
    if args.stage in ['analyse','all']:run_analysis()

if __name__=='__main__':main()
