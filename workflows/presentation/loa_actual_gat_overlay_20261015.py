"""Actual legacy-GAT probability/class overlays, with paired completeness cuts."""
import csv
import json
from pathlib import Path
import numpy as np
from loa_tng_abacus_overlay_20261015 import HOME,KEYS,PROPLABELS,property_data,raw_statistics,sha,dump
from legacy_selection_tests_20261014 import completeness,apply_limit
from desi_density_properties_20261014 import estimate

TRAIN=HOME/'legacy_gat_density_20261015'
OUT=HOME/'loa_actual_gat_overlay_20261015'
LABELS={'abacus':'Abacus P12-A · R7','tng_r2_dense':'GAT · R2 · normal',
        'tng_r2_thin':'GAT · R2 · thinned','tng_r7_dense':'GAT · R7 · normal','tng_r7_thin':'GAT · R7 · thinned'}
WINDOWS=[('full',.2,.3),('lower',.2,.25),('upper',.25,.3)]

def controlled(y,w,c,region):
    try:
        r,_=estimate(y,w,c,region,np.zeros(len(y)),minimum=20,nboot=128)
    except (ValueError,AssertionError) as e:
        return {'unavailable':str(e) or 'Insufficient valid spatial bootstrap draws'}
    r.pop('paired_density_minus_soft_contrast',None)
    error=0.
    for k,v in w.items():
        mean=np.zeros(4)
        for cell,ref in zip(r['common_cell_ids'],r['reference_weights']):
            q=c==cell;mean+=ref*(v[q].T@y[q])/v[q].sum(0)
        error=max(error,float(np.max(abs(mean-r['methods'][k]['mean']))))
    assert error<1e-10
    r['independent_arithmetic_max_error']=error
    return r

def main():
    OUT.mkdir(exist_ok=True)
    training_validation=json.loads((TRAIN/'VALIDATION.json').read_text())
    assert training_validation['verified']
    training_manifest={v['name']:v['sha256'] for v in json.loads((TRAIN/'MANIFEST.json').read_text())['files']}
    source=HOME/'loa_tng_abacus_overlay_20261015/loa_predictions.npz';a=np.load(source)
    assert sha(source)==json.loads((source.parent/'inference_receipt.json').read_text())['prediction_sha256']
    photpath=HOME/'legacy_selection_tests_20261014/p12a_targeting.npz';phot=np.load(photpath)
    assert np.array_equal(a['TARGETID'],phot['TARGETID'])
    good,yy=property_data(a);bright=phot['r']<phot['rlimit'];base=good&bright
    parent=a['CIGALE_MATCHED']&(a['MASS_CG_15']>0)&np.isfinite(yy['mass'])&bright
    curve=completeness(a['z'][parent],yy['mass'][parent],phot['r'][parent],phot['rlimit'][parent],np.linspace(.2,.3,11))
    cuts={'uncut':np.ones(len(base),bool),'complete90':apply_limit(a['z'],yy['mass'],curve)}
    result={'scenarios':{},'paired_cuts':{},'curve':curve,'inputs':{str(source):sha(source),str(photpath):sha(photpath)},
      'training':{},'verification':{},'boundary_sensitivity':{},
      'methods':'Actual four-layer legacy GAT classifiers; one fit per density/radius, probabilities are class scores, not calibrated posteriors. Soft class means and argmax classes compared with corresponding Abacus P_WEB estimators.',
      'selection':'Same previously validated Loa query IDs, CG15 properties and parent mask; no observed-tracer thinning or graph rebuilding after outcome cuts.',
      'normalizations':'legacy_targetfit fits BoxCox to all observed graph nodes and recenters, matching the historical default. frozen applies only the unchanged train-fitted transform. Neither uses observed property labels.',
      'limits':'One seed and 1200 minibatch updates, validation-selected checkpoint; no convergence or DESI calibration claim. Four-hop context can extend beyond the original local-feature boundary mask; a conservative interior sensitivity is recorded.'}
    probs={norm:{'abacus':a['P_WEB'].astype(float)} for norm in ['legacy_targetfit','frozen']}
    for radius in [2,7]:
        for arm in ['dense','thin']:
            name=f'R{radius}_{arm}';p=TRAIN/f'{name}_loa.npz';d=np.load(p);k=f'tng_r{radius}_{arm}'
            assert sha(p)==training_manifest[p.name]
            assert np.array_equal(a['TARGETID'],d['TARGETID'])
            result['inputs'][str(p)]=sha(p)
            result['training'][k]=json.loads((TRAIN/f'{name}_training.json').read_text())
            for norm in probs:
                v=d[norm].astype(float);assert np.isfinite(v).all() and (v>=0).all() and np.allclose(v.sum(1),1,atol=2e-6)
                probs[norm][k]=v
                result['verification'][k+'_'+norm]=float(d[norm+'_layerwise_parity'])
    for norm,pp in probs.items():
        for window,lo,hi in WINDOWS:
            for cut,keep in cuts.items():
                rows=np.flatnonzero(base&keep&(a['z']>=lo)&(a['z']<hi))
                zbin=np.floor((a['z'][rows]-.2+1e-10)/.01).astype(int)
                mb=np.floor((yy['mass'][rows]-9+1e-10)/.1).astype(int)
                _,region=np.unique(a['region'][rows],return_inverse=True);nr=region.max()+1
                draws=np.random.default_rng(20261015).multinomial(nr,np.ones(nr)/nr,size=128)
                soft={k:v[rows] for k,v in pp.items()};hard={k:np.eye(4)[v[rows].argmax(1)] for k,v in pp.items()}
                suites={s:{} for s in ['soft_raw','soft_controlled','hard_raw','hard_controlled']}
                for prop,value in yy.items():
                    y=value[rows].astype(float);cells=zbin if prop=='mass' else zbin+11*mb
                    for estimator,w in [('soft',soft),('hard',hard)]:
                        suites[estimator+'_raw'][prop]=raw_statistics(y,w,region,draws)
                        suites[estimator+'_controlled'][prop]=controlled(y,w,cells,region)
                        # Separate raw reference reconstruction.
                        for k,v in w.items():
                            expected=np.divide(v.T@y,v.sum(0),out=np.full(4,np.nan),where=v.sum(0)>0)
                            assert np.allclose(expected,np.asarray(suites[estimator+'_raw'][prop]['methods'][k]['mean'],float),equal_nan=True,atol=1e-10)
                scenario={'rows':len(rows),'z':[lo,hi],'cut':cut,'suites':suites,
                    'soft_fractions':{k:v.mean(0).tolist() for k,v in soft.items()},'hard_fractions':{k:v.mean(0).tolist() for k,v in hard.items()}}
                key=norm+'_'+window+'_'+cut;result['scenarios'][key]=scenario
                print(key,len(rows),{k:round(v['contrast']*100,2) for k,v in suites['soft_controlled']['quenched'].get('methods',{}).items()},flush=True)
            rows=np.flatnonzero(base&(a['z']>=lo)&(a['z']<hi));_,region=np.unique(a['region'][rows],return_inverse=True)
            cells=np.floor((a['z'][rows]-.2+1e-10)/.01).astype(int)+11*np.floor((yy['mass'][rows]-9+1e-10)/.1).astype(int)
            w={k+'_parent':v[rows] for k,v in pp.items()}
            w.update({k+'_cut':v[rows]*cuts['complete90'][rows,None] for k,v in pp.items()})
            result['paired_cuts'][norm+'_'+window]=controlled(yy['quenched'][rows],w,cells,region)
        # Edges <=2alpha: four message hops plus two hops for neighbor tetra
        # density gives <=12alpha. A conservative local mask-margin diagnostic.
        interior=base&(a['boundary_lower_Mpc']>12*6.906260417222935)
        rows=np.flatnonzero(interior);_,region=np.unique(a['region'][rows],return_inverse=True)
        cells=np.floor((a['z'][rows]-.2+1e-10)/.01).astype(int)+11*np.floor((yy['mass'][rows]-9+1e-10)/.1).astype(int)
        result['boundary_sensitivity'][norm]={'selected_rows':len(rows),'controlled':controlled(yy['quenched'][rows],{k:v[rows] for k,v in pp.items()},cells,region)}
    result['matched_normalizations']={}
    for cut,keep in cuts.items():
        rows=np.flatnonzero(base&keep);_,region=np.unique(a['region'][rows],return_inverse=True)
        cells=np.floor((a['z'][rows]-.2+1e-10)/.01).astype(int)+11*np.floor((yy['mass'][rows]-9+1e-10)/.1).astype(int)
        w={'abacus':probs['frozen']['abacus'][rows]}
        for norm,pp in probs.items():
            w.update({norm+'_'+k:v[rows] for k,v in pp.items() if k!='abacus'})
        result['matched_normalizations'][cut]=controlled(yy['quenched'][rows],w,cells,region)
    result['script_sha256']=sha(__file__)
    result['training_validation_sha256']=sha(TRAIN/'VALIDATION.json')
    dump(OUT/'RESULTS.json',result)
    with (OUT/'curves.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['scenario','suite','property','model','class','mean','lower16','upper84','rows'])
        for key,s in result['scenarios'].items():
            for suite,panels in s['suites'].items():
                for prop,r in panels.items():
                    for model,d in r.get('methods',{}).items():
                        for j in range(4):writer.writerow([key,suite,prop,model,j,d['mean'][j],d['interval'][0][j],d['interval'][1][j],r['rows']])
    render(result)

def render(r):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from shared.plot_style import apply_style,ACCENT_COLORS,TEXT_COLOR
    apply_style();colors=[TEXT_COLOR,ACCENT_COLORS['blue'],ACCENT_COLORS['blue'],ACCENT_COLORS['magenta'],ACCENT_COLORS['magenta']]
    with PdfPages(OUT/'actual_GAT_Abacus_overlays.pdf') as pdf:
        for norm in ['legacy_targetfit','frozen']:
            for suite in ['soft_raw','soft_controlled','hard_raw','hard_controlled']:
                for prop in (['quenched'] if suite.startswith('hard') else list(PROPLABELS)):
                    f,axes=plt.subplots(3,2,figsize=(14,13),sharex=True,sharey=True)
                    for row,(window,lo,hi) in enumerate(WINDOWS):
                        for col,cut in enumerate(['uncut','complete90']):
                            ax=axes[row,col];s=r['scenarios'][norm+'_'+window+'_'+cut];p=s['suites'][suite][prop]
                            for j,k in enumerate(KEYS):
                                if 'methods' not in p:continue
                                d=p['methods'][k];x=np.arange(4)+(j-2)*.025;ci=np.asarray(d['interval'],float)
                                ax.plot(x,np.asarray(d['mean'],float),color=colors[j],ls='--' if k.endswith('dense') else '-',marker='o',ms=4,label=LABELS[k])
                                ax.vlines(x,ci[0],ci[1],color=colors[j],alpha=.65)
                            ax.set_title(f'{lo:.2f} ≤ z < {hi:.2f} · '+('No completeness cut' if cut=='uncut' else '90% limiting-mass cut')+f'\nN = {p.get("rows",s["rows"]):,}',fontsize=12)
                            ax.set_xticks(range(4),['Void','Wall','Filament','Cluster']);ax.grid(alpha=.15)
                            if col==0:ax.set_ylabel(PROPLABELS[prop])
                            if 'methods' not in p:ax.text(.5,.5,'No shared four-class support',transform=ax.transAxes,ha='center')
                    f.suptitle(f'Actual TNG GAT versus Abacus P12-A · {prop}\n{suite} · '+('Historical observed-feature normalization' if norm=='legacy_targetfit' else 'Frozen training-feature normalization'),fontsize=17)
                    handles,labels=axes[0,0].get_legend_handles_labels()
                    if not handles:
                        from matplotlib.lines import Line2D
                        handles=[Line2D([],[],color=colors[j],ls='--' if k.endswith('dense') else '-') for j,k in enumerate(KEYS)];labels=[LABELS[k] for k in KEYS]
                    f.legend(handles,labels,loc='lower center',ncol=3,bbox_to_anchor=(.5,.035))
                    f.text(.5,.012,'One training seed per GAT; 128 sky-block resamples (16–84%). Class scores are not calibrated observational truth.',ha='center',fontsize=10)
                    f.subplots_adjust(top=.90,bottom=.12,hspace=.30,wspace=.12);pdf.savefig(f)
                    if prop=='quenched':f.savefig(OUT/f'{norm}_{suite}_quenched.png',dpi=130)
                    plt.close(f)

if __name__=='__main__':main()
