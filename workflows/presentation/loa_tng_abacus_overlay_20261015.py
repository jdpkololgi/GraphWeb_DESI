"""Apply four frozen TNG model families to the same Loa sample as saved P12-A.

Bounded, single-thread observational comparison. No model fitting, survey
thinning, density rescaling, or production VAC modification.
"""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
HOME=Path('/global/u2/d/dkololgi/fmsu_tokyo_20261010')
SOURCE=HOME/'raw_environment_checks_20261014'
TRAIN=HOME/'tng_density_control_20261015'
OUT=HOME/'loa_tng_abacus_overlay_20261015'
VAC=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1/DESI_LOA_P12A_HALO48_FULL_SURVEY_VAC.fits')
KEYS=['abacus','tng_r2_dense','tng_r2_thin','tng_r7_dense','tng_r7_thin']
LABELS={'abacus':'Abacus P12-A · R7','tng_r2_dense':'TNG · R2 · normal','tng_r2_thin':'TNG · R2 · thinned','tng_r7_dense':'TNG · R7 · normal','tng_r7_thin':'TNG · R7 · thinned'}
PROPLABELS={'mass':r'Mean log$_{10}$ stellar mass [M$_\odot$]',
 'sfr':r'Mean log$_{10}$ SFR [M$_\odot$ yr$^{-1}$]',
 'ssfr':r'Mean log$_{10}$ sSFR [yr$^{-1}$]',
 'quenched':'Quenched fraction', 'age':r'Mean log$_{10}$ fitted age [Gyr]',
 'av':r'Mean dust attenuation $A_V$ [mag]', 'fagn':'Mean fitted AGN fraction'}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def dump(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def unit(ra,dec):
    a=np.deg2rad(ra);b=np.deg2rad(dec)
    return np.column_stack([np.cos(a)*np.cos(b),np.sin(a)*np.cos(b),np.sin(b)])

def prepare(out):
    import fitsio
    from scipy.spatial import cKDTree
    from astropy.cosmology import Planck18
    import torch
    from torch import nn
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    provenance=json.loads((SOURCE/'measurement_provenance.json').read_text())
    assert sha(SOURCE/'raw_environment_catalogue.npz')==provenance['catalogue_sha256']
    a=np.load(SOURCE/'raw_environment_catalogue.npz')
    c=json.loads((TRAIN/'tng_density_control_20261015.json').read_text())
    manifest=json.loads((TRAIN/'MANIFEST.json').read_text())
    expected={v['name']:v['sha256'] for v in manifest['files']}
    cols=['TARGETID','RA','DEC','Z']
    d=fitsio.read(VAC,columns=cols)
    assert len(d)==provenance['source_receipt']['rows']
    qrows=a['catalogue_rows'];q=d[qrows]
    for key,qkey in [('TARGETID','TARGETID'),('RA','RA'),('DEC','DEC'),('Z','z')]:
        assert np.array_equal(q[key],a[qkey]),key
    assert len(np.unique(q['TARGETID']))==len(q)
    h=float(Planck18.h);zg=np.linspace(.1,.6,10001);rg=Planck18.comoving_distance(zg).value*h
    radius=np.interp(d['Z'],zg,rg);xyz=unit(d['RA'],d['DEC'])*radius[:,None];qx=xyz[qrows]
    pad=64.;lo=qx.min(0)-pad;hi=qx.max(0)+pad
    keep=np.all((xyz>=lo)&(xyz<=hi),axis=1);tracer_rows=np.flatnonzero(keep)
    tx=xyz[keep];ids=d['TARGETID'][keep];qr=radius[qrows]
    assert len(np.unique(ids))==len(ids)
    tree=cKDTree(tx)
    counts=np.stack([tree.query_ball_point(qx,r,return_length=True,workers=1)-1 for r in c['feature_apertures_hmpc']],axis=1)
    distance,neighbours=tree.query(qx,k=41,workers=1)
    # Exclude the actual query TARGETID, rather than assuming it wins a distance tie.
    isself=ids[neighbours]==q['TARGETID'][:,None]
    assert (isself.sum(1)==1).all()
    other=np.where(isself,np.inf,distance)
    other.sort(axis=1);nearest=other[:,np.asarray(c['feature_neighbour_ranks'])-1]
    assert nearest[:,-1].max()<pad
    assert np.array_equal(counts[:,:4],a['counts'])
    brute_ids=np.linspace(0,len(q)-1,12).astype(int)
    for i in brute_ids:
        dd=np.linalg.norm(tx-qx[i],axis=1);dd[ids==q['TARGETID'][i]]=np.inf
        assert np.array_equal([(dd<=r).sum() for r in c['feature_apertures_hmpc']],counts[i])
        assert np.allclose(np.sort(dd)[np.asarray(c['feature_neighbour_ranks'])-1],nearest[i],atol=1e-10,rtol=0)
    features=np.column_stack([np.log1p(counts),np.log(np.maximum(nearest,1e-6))]).astype(np.float32)
    assert np.isfinite(features).all()
    radedge=np.minimum(qr-np.interp(.15,zg,rg),np.interp(.55,zg,rg)-qr)
    needed=np.maximum(14.,nearest[:,-1])
    geometric=a['valid_sphere']&(a['boundary_lower_Mpc']*h>needed)&(radedge>needed)
    output={k:a[k] for k in a.files}
    output.update(features=features,feature_counts=counts,feature_nearest=nearest,geometric_complete=geometric)
    receipt={'source_cache_sha256':sha(SOURCE/'raw_environment_catalogue.npz'),'source_provenance':provenance,
      'training_manifest_sha256':sha(TRAIN/'MANIFEST.json'),'query_rows':len(q),'tracer_rows':len(tx),
      'tracer_fields_sha256':hashlib.sha256(d[tracer_rows].tobytes()).hexdigest(),
      'coordinates':f'Planck18 comoving h^-1 Mpc, h={h}; redshift space, no periodic wrapping',
      'feature_contract':c['feature_apertures_hmpc']+c['feature_neighbour_ranks'],
      'geometric_complete_rows':int(geometric.sum()),'max_d40_hmpc':float(nearest[:,-1].max()),
      'brute_force_queries':brute_ids.tolist(),'first_four_counts_equal_cached':True,
      'checkpoints':{},'feature_shift':{},'observed_tracers_thinned':False,
      'note':'Normal/thinned identifies training populations. All models see the same observed Loa features with their own unchanged training scalers.'}
    for r in [2,7]:
        for arm in ['dense','thin']:
            family=f'tng_r{r}_{arm}';pred=[]
            for seed in c['seeds']:
                name=f'R{r}_dense_seed{seed}' if arm=='dense' else f'R{r}_thin_{seed}_seed{seed}'
                path=TRAIN/(name+'.pt');digest=sha(path);assert digest==expected[path.name]
                ck=torch.load(path,map_location='cpu',weights_only=False)
                x=(features-ck['x_mean'])/ck['x_std']
                model=nn.Sequential(nn.Linear(11,64),nn.ReLU(),nn.Linear(64,64),nn.ReLU(),nn.Linear(64,3))
                model.load_state_dict(ck['state']);model.eval()
                with torch.no_grad():
                    z=np.concatenate([model(torch.from_numpy(v)).numpy() for v in np.array_split(x,8)])*ck['y_std']+ck['y_mean']
                eig=np.column_stack([z[:,0],z[:,0,None]+np.cumsum(np.logaddexp(0,z[:,1:]),axis=1)]).astype(np.float32)
                assert np.isfinite(eig).all() and (np.diff(eig,axis=1)>=0).all()
                output[name]=eig;pred.append(eig);receipt['checkpoints'][name]={'path':str(path),'sha256':digest}
                receipt['feature_shift'][name]={'any_feature_outside_5_training_sd_fraction':float((np.abs(x[geometric])>5).any(1).mean()),'median_standardized_features':np.median(x[geometric],axis=0).tolist()}
            output[family]=np.mean(pred,axis=0)
            print('Inferred',family,len(q),flush=True)
    output['abacus']=a['EIGENVALUE_MEAN'].astype(np.float32)
    np.savez_compressed(out/'loa_predictions.npz',**output)
    receipt['prediction_sha256']=sha(out/'loa_predictions.npz')
    receipt['script_sha256']=sha(__file__)
    dump(out/'inference_receipt.json',receipt)
    print('Features and inference complete; geometric queries',geometric.sum(),flush=True)

def memberships(values,strata):
    from desi_raw_environment_plots_20261014 import tied_quartiles
    return tied_quartiles(values,strata)[0]

def property_data(a):
    lm=np.log10(np.maximum(a['MASS_CG_15'],1e-30));ls=np.log10(np.maximum(a['SFR_CG_15'],1e-30))
    good=a['geometric_complete']&a['CIGALE_MATCHED']&(np.char.strip(a['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')
    for name in ['MASS_CG_15','SFR_CG_15','AGE_CG_15']:good&=np.isfinite(a[name])&(a[name]>0)
    good&=np.isfinite(a['AV_CG_15'])&(a['AV_CG_15']>=0)&np.isfinite(a['FAGN_CG_15'])&(a['FAGN_CG_15']>=0)&(a['FAGN_CG_15']<=1)&(lm>=9)&(lm<12)
    y={'mass':lm,'sfr':ls,'ssfr':ls-lm,'quenched':(ls-lm<-11).astype(float),
       'age':np.log10(np.maximum(a['AGE_CG_15'],1e-30)/1000),'av':a['AV_CG_15'],'fagn':a['FAGN_CG_15']}
    return good,y

def raw_statistics(y,weights,region,draws):
    nr=draws.shape[1];out={}
    for name,w in weights.items():
        count=np.stack([np.bincount(region,weights=w[:,j],minlength=nr) for j in range(4)])
        sums=np.stack([np.bincount(region,weights=w[:,j]*y,minlength=nr) for j in range(4)])
        totals=count.sum(1);means=np.divide(sums.sum(1),totals,out=np.full(4,np.nan),where=totals>0)
        cb=draws@count.T;sb=draws@sums.T
        bs=np.divide(sb,cb,out=np.full_like(sb,np.nan),where=cb>0)
        intervals=[]
        for j in range(4):
            vals=bs[:,j];vals=vals[np.isfinite(vals)]
            intervals.append(np.quantile(vals,[.16,.84]).tolist() if len(vals)>=.8*len(draws) else [None,None])
        out[name]={'mean':[float(x) if np.isfinite(x) else None for x in means],
                   'interval':np.asarray(intervals,dtype=object).T.tolist(),
                   'group_weight':totals.tolist(),
                   'contrast':float(means[3]-means[0]) if np.isfinite(means[[0,3]]).all() else None}
    return {'rows':len(y),'methods':out,'sky_blocks':nr,'bootstrap_valid':len(draws)}

def analyse(out):
    from desi_density_properties_20261014 import estimate
    a=np.load(out/'loa_predictions.npz');receipt=json.loads((out/'inference_receipt.json').read_text())
    assert sha(out/'loa_predictions.npz')==receipt['prediction_sha256']
    good,yall=property_data(a);ix=np.flatnonzero(good)
    y={k:v[ix].astype(float) for k,v in yall.items()}
    zb=np.floor((a['z'][ix]-.2+1e-10)/.01).astype(int)
    mb=np.floor((y['mass']-9+1e-10)/.1).astype(int)
    _,regions=np.unique(a['region'][ix],return_inverse=True);nr=regions.max()+1
    rng=np.random.default_rng(20261015);draws=rng.multinomial(nr,np.ones(nr)/nr,size=256)
    quart={k:memberships(a[k][ix].sum(1),zb) for k in KEYS}
    classes={k:np.eye(4)[(a[k][ix]>=.2).sum(1)] for k in KEYS}
    suites={'rank_raw':{},'rank_controlled':{},'class_raw':{},'class_soft_abacus':{}}
    softclasses=dict(classes);softclasses['abacus']=a['P_WEB'][ix].astype(float)
    assert np.allclose(softclasses['abacus'].sum(1),1,atol=1e-6)
    support={};validation={};seed_stats={}
    for name,values in y.items():
        suites['rank_raw'][name]=raw_statistics(values,quart,regions,draws)
        suites['class_raw'][name]=raw_statistics(values,classes,regions,draws)
        suites['class_soft_abacus'][name]=raw_statistics(values,softclasses,regions,draws)
        cells=zb if name=='mass' else zb+11*mb
        r,common=estimate(values,quart,cells,regions,np.zeros(len(ix)),minimum=20,nboot=256)
        r.pop('paired_density_minus_soft_contrast',None)
        suites['rank_controlled'][name]=r;support[name]=ix[common]
        # Independent direct row-wise reconstruction at the saved common mix.
        maxerr=0.
        for model,w in quart.items():
            means=np.zeros(4)
            for cell,ref in zip(r['common_cell_ids'],r['reference_weights']):
                ii=cells==cell
                means+=ref*np.sum(w[ii]*values[ii,None],axis=0)/w[ii].sum(0)
            maxerr=max(maxerr,float(np.max(abs(means-r['methods'][model]['mean']))))
        assert maxerr<1e-10;validation[name]=maxerr
    # Class control cannot invent absent classes; record common support explicitly.
    try:
        class_controlled={}
        for name,values in y.items():
            cells=zb if name=='mass' else zb+11*mb
            r,common=estimate(values,classes,cells,regions,np.zeros(len(ix)),minimum=20,nboot=256)
            r.pop('paired_density_minus_soft_contrast',None);class_controlled[name]=r
        suites['class_controlled']=class_controlled
        class_control_status='available'
    except ValueError:
        class_control_status='No shared four-class support across all five models; raw class curves shown with gaps for empty groups.'
    # Seed sensitivity on the same retained primary support, without refitting.
    for r in [2,7]:
        for arm in ['dense','thin']:
            family=f'tng_r{r}_{arm}';vals=[]
            for seed in [20261015,20261016,20261017]:
                ck=f'R{r}_dense_seed{seed}' if arm=='dense' else f'R{r}_thin_{seed}_seed{seed}'
                w=memberships(a[ck][ix].sum(1),zb);ref=suites['rank_controlled']['quenched'];means=np.zeros(4);cells=zb+11*mb
                for cell,mix in zip(ref['common_cell_ids'],ref['reference_weights']):
                    ii=cells==cell;den=w[ii].sum(0)
                    if (den==0).any():raise ValueError('Seed lacks ensemble common support')
                    means+=mix*(w[ii].T@y['quenched'][ii])/den
                vals.append(float(means[-1]-means[0]))
            seed_stats[family]={'quenched_quartile_contrasts':vals,'minmax':[min(vals),max(vals)]}
    fractions={k:classes[k].mean(0).tolist() for k in KEYS}
    # Saved Abacus posterior probabilities are not labels of posterior means.
    soft=raw_statistics(y['quenched'],{'abacus_probabilities':a['P_WEB'][ix].astype(float)},regions,draws)
    result={'suites':suites,'model_labels':LABELS,'property_labels':PROPLABELS,'input_property_rows':len(ix),
      'class_fractions':fractions,'class_control_status':class_control_status,'seed_sensitivity':seed_stats,
      'abacus_probability_weighted_quenched_check':soft,'direct_controlled_arithmetic_max_errors':validation,
      'definition':'Five same-galaxy model families. TNG predictions average physical eigenvalues from three separately frozen fits. Abacus uses saved posterior-mean eigenvalues. Ranks are per-model trace quartiles within dz=.01. Class plots count eigenvalues >=.2, including for Abacus.',
      'control':'Identical common .1-dex mass/.01-z cells and pooled reference mix across all five models. Mass panel controls z only. Minimum Kish effective count 20 per group/cell. Raw plots use identical input galaxies with no reweighting.',
      'uncertainty':'256 paired nside8 sky-block resamples, fixed groups/control weights, 16-84 percentile intervals. Conditional sampling variation only; excludes property errors, model transfer, posterior uncertainty and fit-seed uncertainty.',
      'selection':'Same Loa wedge RA120-160 Dec14.5-30.6 z.20-.30; conservative mask/radial margin larger than max(14 h^-1 Mpc,40th-neighbour distance); same finite valid positive mass/SFR/age, AV>=0, FAGN in[0,1], CG15 galaxy match,9<=logM<12.',
      'scientific_status':'Descriptive zero-shot comparison; no observed environment truth or DESI calibration claim.'}
    dump(out/'RESULTS.json',result)
    np.savez_compressed(out/'analysis_support.npz',property_ids=ix,**support)
    with (out/'curves.csv').open('w') as file:
        writer=csv.writer(file);writer.writerow(['suite','property','model','group','mean','lower16','upper84','rows'])
        for suite,panels in suites.items():
            for prop,r in panels.items():
                for model,d in r['methods'].items():
                    for j in range(4):writer.writerow([suite,prop,model,j+1,d['mean'][j],d['interval'][0][j],d['interval'][1][j],r['rows']])
    print(json.dumps({'rows':len(ix),'classes':fractions,'controlled_quenched':suites['rank_controlled']['quenched'],'seed':seed_stats},indent=2),flush=True)

def render(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from shared.plot_style import apply_style,finalize_axes,ACCENT_COLORS,TEXT_COLOR,COSMIC_WEB_COLORS,CLASS_ORDER
    from PIL import Image
    apply_style();r=json.loads((out/'RESULTS.json').read_text())
    style={'abacus':(TEXT_COLOR,'D','-'),'tng_r2_dense':(ACCENT_COLORS['blue'],'o','--'),
      'tng_r2_thin':(ACCENT_COLORS['blue'],'s','-'),'tng_r7_dense':(ACCENT_COLORS['magenta'],'^','--'),
      'tng_r7_thin':(ACCENT_COLORS['magenta'],'v','-')}
    pdf=PdfPages(out/'LOA_TNG_Abacus_comparison.pdf')
    def save(fig,name):
        for ext in ['png','pdf','svg']:fig.savefig(out/f'{name}.{ext}',dpi=180,bbox_inches='tight')
        pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
        with Image.open(out/f'{name}.png') as im:
            im.thumbnail((1800,1200));im.save(out/f'{name}_preview.png')
    def curves(suite,compact=False):
        properties=['quenched','ssfr','age'] if compact else list(PROPLABELS)
        f,axes=plt.subplots(1,3,figsize=(16,6.8)) if compact else plt.subplots(2,4,figsize=(21,11))
        f.subplots_adjust(left=.07,right=.98,bottom=.21 if compact else .16,top=.74 if compact else .79,wspace=.36,hspace=.45)
        isclass=suite.startswith('class');controlled=suite.endswith('controlled');softabacus=suite=='class_soft_abacus'
        title='Loa properties versus inferred web class' if isclass else 'Loa properties versus inferred environment rank'
        f.suptitle(title,fontsize=22,y=.98)
        f.text(.5,.925,'Same DESI Loa galaxies · 0.20 ≤ z < 0.30 · CIGALE 15-band properties',ha='center',fontsize=14)
        for ax,prop in zip(np.asarray(axes).flat,properties):
            panel=r['suites'][suite][prop]
            for i,key in enumerate(KEYS):
                d=panel['methods'][key];mu=np.asarray(d['mean'],float);ci=np.asarray(d['interval'],float);color,marker,ls=style[key]
                x=np.arange(4)+(i-2)*.035
                label='Abacus P12-A · posterior weights' if key=='abacus' and softabacus else LABELS[key]
                ax.plot(x,mu,color=color,marker=marker,ls=ls,lw=1.8,markersize=5,label=label)
                ax.vlines(x,ci[0],ci[1],color=color,lw=1,alpha=.8)
            if isclass:ax.set_xticks(range(4),['Void','Wall','Filament','Cluster'])
            else:ax.set_xticks(range(4),['Q1','Q2','Q3','Q4'])
            finalize_axes(ax,title='',xlabel='Inferred web class' if isclass else 'Low → high inferred trace quartile',ylabel=PROPLABELS[prop],legend=False)
            ax.grid(axis='y',alpha=.15)
            ax.text(.025,.97,f"N = {panel['rows']:,}",transform=ax.transAxes,va='top',fontsize=10)
            if prop=='quenched':
                if controlled:
                    bound=np.asarray([v['interval'] for v in panel['methods'].values()],float)
                    ax.set_ylim(max(0,float(np.nanmin(bound))-.025),min(1,float(np.nanmax(bound))+.025))
                else:ax.set_ylim(0,1)
            else:ax.margins(y=.30)
        if not compact:
            last=np.asarray(axes).flat[-1];last.axis('off')
            note='Same property sample for all models.\n\nNormal / thinned refer to TNG\ntraining density, not Loa thinning.\n\nAbacus: saved P12-A predictions.\nTNG: mean of three frozen fits.\n\nR2 / R7 are target smoothing\nscales in comoving Mpc/h.'
            last.text(0,.95,note,va='top',linespacing=1.6,fontsize=12)
        handles,labels=np.asarray(axes).flat[0].get_legend_handles_labels()
        f.legend(handles,labels,ncol=3,loc='upper center',bbox_to_anchor=(.5,.89),fontsize=11)
        control='Common mass + redshift mix; stellar-mass panel controls redshift only.' if controlled else 'Raw selected-sample means; stellar-mass and redshift mixtures are not balanced.'
        f.text(.5,.105 if compact else .085,control,ha='center',fontsize=12)
        foot='Class gaps indicate empty groups; Abacus classes use posterior-mean eigenvalues, not posterior probabilities.' if isclass else 'Each model defines its own trace quartiles within Δz = 0.01; equal rank is not equal physical environment.'
        if softabacus:foot='Abacus: posterior-probability weights. TNG: classes of mean predicted eigenvalues. Empty classes remain gaps.'
        f.text(.5,.065 if compact else .05,foot,ha='center',fontsize=11)
        f.text(.5,.025,'Bars: 16–84% paired sky-block bootstrap · zero-shot transfer; inferred environments are not calibrated observational truth',ha='center',fontsize=10)
        save(f,('summary_' if compact else '')+suite)
    for suite in ['rank_raw','rank_controlled','class_raw','class_soft_abacus']:
        curves(suite,False)
        if suite!='rank_raw':curves(suite,True)
    if 'class_controlled' in r['suites']:curves('class_controlled',False)
    f,ax=plt.subplots(figsize=(13,6));bottom=np.zeros(5)
    for j,cls in enumerate(CLASS_ORDER):
        height=np.array([r['class_fractions'][k][j] for k in KEYS]);ax.bar(np.arange(5),height,bottom=bottom,color=COSMIC_WEB_COLORS[cls],label=cls.title());bottom+=height
    ax.set_xticks(np.arange(5),['Abacus\nP12-A R7','TNG R2\nnormal','TNG R2\nthinned','TNG R7\nnormal','TNG R7\nthinned'])
    finalize_axes(ax,'Class occupancy on the same Loa galaxies','Model / training population','Fraction of selected galaxies')
    ax.set_ylim(0,1);f.subplots_adjust(bottom=.23,top=.85)
    f.text(.5,.045,f"N = {r['input_property_rows']:,} · classes from mean eigenvalues at threshold 0.2 · identical observed tracers",ha='center',fontsize=12)
    save(f,'class_occupancy');pdf.close()

def verify(out):
    import fitsio
    import torch
    torch.set_num_threads(1)
    a=np.load(out/'loa_predictions.npz');r=json.loads((out/'RESULTS.json').read_text());receipt=json.loads((out/'inference_receipt.json').read_text())
    good,yall=property_data(a);ix=np.flatnonzero(good);checks=0;max_replay=0.;max_mean=0.
    assert len(ix)==r['input_property_rows'];checks+=1
    sample=np.linspace(0,len(a['TARGETID'])-1,24).astype(int)
    cols=['TARGETID','Z','EIGENVALUE_MEAN','P_WEB']
    source=fitsio.read(VAC,rows=a['catalogue_rows'][sample],columns=cols)
    for key,other in [('TARGETID','TARGETID'),('Z','z'),('EIGENVALUE_MEAN','abacus'),('P_WEB','P_WEB')]:
        assert np.array_equal(source[key],a[other][sample]);checks+=1
    cg=fitsio.read(receipt['source_provenance']['source_receipt']['path'],rows=a['catalogue_rows'][sample],columns=['TARGETID','MASS_CG_15','SFR_CG_15','AGE_CG_15','AV_CG_15','FAGN_CG_15'])
    for key in cg.dtype.names:
        assert np.array_equal(cg[key],a[key][sample],equal_nan=True);checks+=1
    for name,meta in receipt['checkpoints'].items():
        ck=torch.load(meta['path'],map_location='cpu',weights_only=False);assert sha(meta['path'])==meta['sha256'];checks+=1
        z=((a['features'][sample]-ck['x_mean'])/ck['x_std']).astype(np.float32)
        for layer in [0,2,4]:
            z=z@ck['state'][f'{layer}.weight'].numpy().T+ck['state'][f'{layer}.bias'].numpy()
            if layer<4:z=np.maximum(z,0)
        z=z*ck['y_std']+ck['y_mean']
        eig=np.column_stack([z[:,0],z[:,0,None]+np.cumsum(np.logaddexp(0,z[:,1:]),axis=1)])
        err=float(np.max(abs(eig-a[name][sample])));assert np.allclose(eig,a[name][sample],atol=1e-5,rtol=1e-5),(name,err)
        max_replay=max(max_replay,err);checks+=1
    zb=np.floor((a['z'][ix]-.2+1e-10)/.01).astype(int)
    for suite in ['rank_raw','class_raw','class_soft_abacus']:
        for model in KEYS:
            if suite=='rank_raw':w=memberships(a[model][ix].sum(1),zb)
            elif suite=='class_soft_abacus' and model=='abacus':w=a['P_WEB'][ix].astype(float)
            else:w=np.eye(4)[(a[model][ix]>=.2).sum(1)]
            assert np.allclose(w.sum(1),1,atol=1e-6);checks+=1
            for prop,y in yall.items():
                expected=r['suites'][suite][prop]['methods'][model]['mean']
                for j in range(4):
                    if w[:,j].sum()==0:assert expected[j] is None
                    else:
                        err=float(abs(np.average(y[ix],weights=w[:,j])-expected[j]));assert err<1e-10
                        max_mean=max(max_mean,err)
                    checks+=1
    for model in KEYS:
        assert np.isfinite(a[model]).all() and (np.diff(a[model],axis=1)>=0).all();checks+=1
        assert np.isclose(sum(r['class_fractions'][model]),1);checks+=1
    # Focused null/tie tests for the reused estimator and new raw curves.
    w=memberships(np.zeros(80),np.repeat([0,1],40));assert np.allclose(w,.25);checks+=1
    raw=raw_statistics(np.full(80,3.),{'constant':w},np.arange(80)%4,np.ones((10,4),int))
    assert np.allclose(raw['methods']['constant']['mean'],3);checks+=1
    checks+=len(r['direct_controlled_arithmetic_max_errors'])
    result={'passed':True,'checks':checks,'checkpoint_numpy_replay_max_error':max_replay,'raw_direct_mean_max_error':max_mean,
      'controlled_direct_mean_max_error':max(r['direct_controlled_arithmetic_max_errors'].values()),
      'source_spot_checks':len(sample),'inference_count_checks':12,'prediction_sha256':sha(out/'loa_predictions.npz'),
      'script_sha256':sha(__file__),'estimator_sha256':sha(Path(__file__).with_name('desi_density_properties_20261014.py')),
      'tie_helper_sha256':sha(Path(__file__).with_name('desi_raw_environment_plots_20261014.py'))}
    dump(out/'VALIDATION.json',result);print(json.dumps(result,indent=2),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','analyse','render','verify']);p.add_argument('--out',type=Path,default=OUT);args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True);globals()[args.stage](args.out)

if __name__=='__main__':main()
