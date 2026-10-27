"""Package bounded field displays and property aggregates for private Page embeds."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='4'
import base64
import gzip
import json
from pathlib import Path
import struct
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parent))
from loa_model_pilot_20261015 import ROOT, read, sha, save_json, require_compute
import numpy as np
from scipy.ndimage import map_coordinates
from workflows.sbi import e2e_coupled_operators as op
from workflows.sbi import e2e_footprint_tidal_operator as tidal

HERE=Path(__file__).resolve().parent/'loa_model_page_20261015'
MODELS=[f'{a}_{p}_{s}' for a,p in [('legacy','v'),('unet','x'),('transformer','v')] for s in (17,29)]
NAMES=['Earlier CFM · seed 17','Earlier CFM · seed 29','Expanded U-Net · seed 17','Expanded U-Net · seed 29','Expanded transformer · seed 17','Expanded transformer · seed 29']
MASS=np.arange(9,12.01,.5); Z=np.array([.15,.25,.35,.45,.55]); R=np.array([0,17,18,18.5,19,19.5])
PROPS=['log stellar mass','log SFR','log sSFR','log stellar age','Dust attenuation','AGN fraction','Quenched fraction']
UNITS=['log10(M☉)','log10(M☉ / yr)','log10(1 / yr)','log10(Gyr)','mag','fraction','fraction']

def pack(data):
    return base64.b64encode(gzip.compress(data,compresslevel=9,mtime=0)).decode()

def build():
    require_compute(); start=time.monotonic()
    prep=read(ROOT/'PREPARED.json'); complete=read(ROOT/'INFERENCE_COMPLETE.json')
    assert complete['passed'] and len(complete['records'])==12
    output=ROOT/'page_v3';output.mkdir(exist_ok=False)
    field={'shape':[16,8,8],'cell':13.532,'models':MODELS,'names':NAMES,'caps':[],'quantization':{'log2_offset':8,'log2_scale':4096,'probability_scale':32,'entropy_scale':255,'redshift_scale':65535}}
    fields={}; field_bytes=0; records=[]; checks=[]; populations=[]
    for capindex,case in enumerate(prep['cases']):
        chunks=[]
        cap=case['cap']; directory=ROOT/cap
        with np.load(directory/'conditions.npz') as f:
            support=f['joint'][1,32:64,24:40,24:40]; redshift=f['joint'][11,32:64,24:40,24:40]
        mask=(support.reshape(16,2,8,2,8,2).min((1,3,5))>=.5).astype('u1')
        zv=op.mean_pool(redshift,2)
        with np.load(directory/'positions.npz') as f: position=f['position']; xyz=f['xyz']; ids=f['targetid']
        with np.load(directory/'properties.npz') as f: prop={k:f[k] for k in f.files}
        np.testing.assert_array_equal(ids,prop['TARGETID'])
        gs=map_coordinates(support,position.T,order=1,mode='nearest',prefilter=False)
        # Context display uses an ID-ordered deterministic sample, never sends IDs.
        order=np.argsort(ids); sample=order[np.linspace(0,len(order)-1,min(400,len(order)),dtype=int)]
        points=np.rint(position[sample]*1024).astype('<i2')
        samplez=np.rint(prop['Z'][sample]*65535).astype('<u2')
        chunks.extend([mask.tobytes(),np.rint(zv*65535).astype('<u2').tobytes(),points.tobytes(),samplez.tobytes()])
        meta={'name':cap,'origin':case['core_origin_mpc_h'],'pointCount':len(sample),'galaxies':len(ids),
              'zRange':case['z_range'],'supportedDisplayCells':int(mask.sum()),'models':[]}
        good=(prop['SUPPORTED'].astype(bool)&prop['CIGALE_MATCHED'].astype(bool)&~prop['CIGALE_AMBIGUOUS'].astype(bool)&
              (np.char.strip(prop['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')&(gs>=.5))
        massraw=prop['MASS_CG_15']; sfrraw=prop['SFR_CG_15']; ageraw=prop['AGE_CG_15']
        mass=np.full(len(ids),np.nan); sfr=mass.copy(); age=mass.copy()
        for dest,raw in [(mass,massraw),(sfr,sfrraw),(age,ageraw)]:
            ok=np.isfinite(raw)&(raw>0);dest[ok]=np.log10(raw[ok])
        age-=3
        ssfr=sfr-mass; quenched=np.where(np.isfinite(ssfr),(ssfr<-11).astype(float),np.nan)
        y=np.column_stack([mass,sfr,ssfr,age,prop['AV_CG_15'],prop['FAGN_CG_15'],quenched])
        y[~(np.isfinite(y[:,4])&(y[:,4]>=0)),4]=np.nan
        y[~(np.isfinite(y[:,5])&(y[:,5]>=0)&(y[:,5]<=1)),5]=np.nan
        mb=np.searchsorted(MASS,mass,side='right')-1; zb=np.searchsorted(Z,prop['Z'],side='right')-1
        rb=np.searchsorted(R,prop['r_bright'],side='right')-1
        good &= np.isfinite(mass)&np.isfinite(prop['r_bright'])&(mb>=0)&(mb<6)&(zb>=0)&(zb<4)&(rb>=0)&(rb<5)
        populations.append({'cap':cap,'coreRows':len(ids),'eligibleMassBrightnessCG15':int(good.sum()),
            'propertyValid':{name:int((good&np.isfinite(y[:,j])).sum()) for j,name in enumerate(PROPS)}})
        weights=[]
        for modelid in MODELS:
            result=directory/modelid; receipt=read(result/'COMPLETE.json')
            if sha(result/'draws.npz')!=receipt['output_sha256']:raise ValueError('Draw checksum mismatch')
            with np.load(result/'draws.npz') as f:
                rho=f['rho']; tensors=f['tensor']; ge=f['galaxy_eigen']
            pooled=np.stack([op.mean_pool(v.astype(float),2) for v in rho])
            logs=np.log2(pooled)
            if not np.isfinite(logs).all() or logs.min() < -8 or logs.max() >= 8:
                raise ValueError('Display quantization range would clip density')
            q=np.rint((logs+8)*4096).astype('<u2')
            reconstructed=2**(q.astype(float)/4096-8)
            rel=float(np.max(abs(reconstructed-pooled)/pooled))
            if rel > .00009:raise ValueError('Density quantization tolerance')
            chunks.append(q.tobytes())
            eig=tidal.eigensystem(tensors)
            for threshold in (0.,.2):
                labels=(eig>threshold).sum(-1)
                p=np.stack([(labels==i).mean(0) for i in range(4)])
                pooledp=np.stack([op.mean_pool(v,2) for v in p])
                pq=np.rint(pooledp*32).astype('u1')
                assert np.array_equal(pq.sum(0),np.full((16,8,8),32))
                np.testing.assert_allclose(pq/32,pooledp,atol=0,rtol=0)
                ent=-np.sum(np.where(p>0,p*np.log(np.maximum(p,1e-30)),0),axis=0)/np.log(4)
                eq=np.rint(op.mean_pool(ent,2)*255).astype('u1')
                chunks.extend([pq.tobytes(),eq.tobytes()])
            meta['models'].append({'id':modelid,'densityRange':[float(pooled.min()),float(pooled.max())],
                'maxDisplayRelativeError':rel,'draws':4})
            labels=(ge>.2).sum(-1); weights.append(np.stack([(labels==i).mean(0) for i in range(4)],axis=1))
        weights.append(prop['P_WEB'])
        # Aggregate only: no observed identifiers/individual CIGALE values shipped.
        cell=mb*20+zb*5+rb
        for mi,w in enumerate(weights):
            np.testing.assert_allclose(w[good].sum(1),1,rtol=0,atol=1e-5)
            for pi in range(7):
                valid=good&np.isfinite(y[:,pi]); expected=int(valid.sum()); total=0
                for env in range(4):
                    ww=w[valid,env].astype(float); cc=cell[valid]; yy=y[valid,pi]
                    n=np.bincount(cc,weights=ww,minlength=120)
                    n2=np.bincount(cc,weights=ww**2,minlength=120)
                    sy=np.bincount(cc,weights=ww*yy,minlength=120)
                    for ci in np.flatnonzero(n>1e-8):
                        records.append((capindex,mi,int(ci//20),int(ci%20//5),int(ci%5),env,pi,float(n[ci]),float(n2[ci]),float(sy[ci])))
                    total+=n.sum()
                assert abs(total-expected)<max(.01,expected*2e-6)
                checks.append({'cap':cap,'model':mi,'property':pi,'rows':expected,'totalWeight':float(total)})
        raw=b''.join(chunks)
        fields['field_'+cap]={**field,'caps':[meta],'data':pack(raw),'bytes':len(raw)}
        field_bytes+=len(raw)
    properties={'models':MODELS+['p12a'],'names':NAMES+['P12-A galaxy posterior'],'properties':PROPS,'units':UNITS,
        'massEdges':MASS.tolist(),'zEdges':Z.tolist(),'rEdges':R.tolist(),'records':len(records),'recordBytes':19,
        'data':pack(b''.join(struct.pack('<7B3f',*r) for r in records)), 'populations':populations,
        'threshold':.2,'selection':'unique unambiguous CIGALE GALAXY, supported native core and support ≥0.5, 9≤logM<12, north-corrected r<19.5; validity is property-specific',
        'uncertainty':'Four field draws; probability-weighted descriptive means, no inferential confidence interval. Effective counts do not account for spatial dependence.'}
    products={**fields,'properties':properties}
    for name,data in products.items():
        template=(HERE/(('field' if name.startswith('field_') else name)+'.html')).read_text()
        html=template.replace('__DATA__',json.dumps(data,separators=(',',':'),ensure_ascii=False,allow_nan=False))
        if len(html.encode())>256*1024:raise ValueError(name+' exceeds Page embed size')
        (output/(name+'.html')).write_text(html)
    save_json(output/'PACKAGING.json',{'passed':True,'source':sha(__file__),'fieldBytes':field_bytes,'propertyRecords':len(records),
        'files':{n: {'bytes':(output/(n+'.html')).stat().st_size,'sha256':sha(output/(n+'.html'))} for n in products},
        'populations':populations,'aggregate_checks':checks,'seconds':time.monotonic()-start,
        'display':'2x spatial pooling; density log quantization error <0.009%; environment expected native-cell fractions exact at 1/32; entropy quantization <=1/510',
        'context_markers':'400 deterministic ID-ordered unclassified context galaxies per cap; IDs omitted'})
    save_json(output/'REFERENCE_ROWS.json',{'firstRows':[list(r) for r in records[:20]],'lastRows':[list(r) for r in records[-20:]]})
    print(json.dumps({'stage':'page_packaged','files':{n:(output/(n+'.html')).stat().st_size for n in products},'records':len(records)}),flush=True)

if __name__=='__main__':
    build()
    import subprocess
    subprocess.run([sys.executable, str(Path(__file__).with_name('test_loa_model_page_20261015.py'))], check=True)
