"""Render audited environment-resolved colour–mass diagrams and trend products."""
import os
os.environ['MPLBACKEND']='Agg'
import argparse
import csv
import json
from pathlib import Path
import socket
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from run import OUT, WINDOWS, ENV, masks, sha, write

LABELS={'quenched':'Quenched fraction (sSFR < 10⁻¹¹ yr⁻¹)','logmass':'Mean log₁₀ stellar mass [M☉]',
        'logssfr':'Mean log₁₀ sSFR [yr⁻¹; floor −14]','logsfr':'Mean log₁₀ SFR [M☉/yr; floor −14]',
        'colour':'Mean colour [mag]','logage':'Mean log₁₀ fitted age [Gyr]','av':'Fitted dust Aᵥ [mag]','fagn':'Fitted AGN contribution'}
FIG=OUT/'presentation'
XE=np.linspace(10.3,12.5,177);YE=np.linspace(-.5,2.5,181)


def histogram(m,c,w):
    valid=np.isfinite(m)&np.isfinite(c)&np.isfinite(w)&(w>=0)
    shown=valid&(m>=XE[0])&(m<=XE[-1])&(c>=YE[0])&(c<=YE[-1])
    h=np.histogram2d(m[valid],c[valid],bins=(XE,YE),weights=w[valid])[0]
    assert np.isclose(h.sum(),w[shown].sum(),rtol=1e-12,atol=1e-7)
    xi=np.minimum(np.searchsorted(XE,m[shown],side='right')-1,len(XE)-2)
    yi=np.minimum(np.searchsorted(YE,c[shown],side='right')-1,len(YE)-2)
    check=np.bincount(xi*(len(YE)-1)+yi,weights=w[shown],minlength=h.size).reshape(h.shape)
    assert np.allclose(h,check,rtol=1e-12,atol=1e-8)
    return h,dict(weight_total=float(w[valid].sum()),weight_displayed=float(w[shown].sum()),weight_outside=float(w[valid&~shown].sum()))


def draw(ax,h,limit,title,note,norm):
    ax.axvspan(XE[0],limit,color='#e7e7e7',zorder=0)
    mesh=ax.pcolormesh(XE,YE,np.ma.masked_equal(h.T,0),cmap='magma',norm=norm,rasterized=True)
    ax.axvline(limit,color='#333333',ls='--',lw=1)
    ax.axhline(.75,color='#158ac2',ls='--',lw=1)
    ax.axhline(1.5,color='#b43d4b',ls=':',lw=1)
    ax.set(xlim=(XE[0],XE[-1]),ylim=(YE[0],YE[-1]),title=title)
    ax.text(.03,.97,note,transform=ax.transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    return mesh


def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(FIG/f'{name}.{ext}',dpi=170)
    plt.close(fig)


def diagrams(windows):
    import fitsio
    data={};hist={};audit={}
    for key,*_ in WINDOWS:
        for variant in ['observed_original','rest_primary','paired']:
            data[key,variant]=fitsio.read(OUT/f'{variant}_{key}.fits',columns=['MASS_CG_15','GR_OBS','GR_REST01','P_WEB'])
        for frame,variant,col in [('observed','observed_original','GR_OBS'),('rest','rest_primary','GR_REST01')]:
            a=data[key,variant];m=np.log10(a['MASS_CG_15']);hard=a['P_WEB'].argmax(1)
            for method in ['hard','soft']:
                for j in range(4):
                    weights=(hard==j).astype(float) if method=='hard' else a['P_WEB'][:,j].astype(float)
                    h,info=histogram(m,a[col],weights);tag=f'{frame}_{method}_{key}_{j}';hist[tag]=h;audit[tag]=info
        for frame,col in [('observed','GR_OBS'),('rest','GR_REST01')]:
            a=data[key,'paired'];tag=f'paired_{frame}_{key}'
            hist[tag],audit[tag]=histogram(np.log10(a['MASS_CG_15']),a[col],np.ones(len(a)))
    # Same colour scale across observed/rest, all windows, and all classes.
    vmax=max(h.max() for k,h in hist.items() if not k.startswith('paired'))
    norm=LogNorm(vmin=.1,vmax=vmax)
    for frame in ['observed','rest']:
        for method in ['hard','soft']:
            fig,axs=plt.subplots(3,4,figsize=(18,13),sharex=True,sharey=True)
            fig.subplots_adjust(left=.065,right=.91,bottom=.10,top=.90,wspace=.12,hspace=.28)
            for i,(key,*_) in enumerate(WINDOWS):
                w=windows[key];lim=w['old_limit'] if frame=='observed' else w['rest_limit']
                for j,ax in enumerate(axs[i]):
                    tag=f'{frame}_{method}_{key}_{j}';rec=audit[tag]
                    title=f"{ENV[j]} · {w['zlo']:.2f} ≤ z < {w['zhi']:.2f}"
                    note=f"{'N' if method=='hard' else 'Σ p'} = {rec['weight_total']:,.0f}\nlog M★ ≥ {lim:.2f}\nOutside axes: {rec['weight_outside']:.1f}"
                    mesh=draw(ax,hist[tag],lim,title,note,norm)
                    if j==0:ax.set_ylabel('Observed g − r [mag]' if frame=='observed' else 'Rest-frame ⁰·¹(g − r) [mag]')
                    if i==2:ax.set_xlabel('log₁₀ stellar mass / M☉ (CIGALE)')
            fig.colorbar(mesh,cax=fig.add_axes([.93,.18,.012,.64]),label='Count / probability weight per bin (log scale)')
            frame_label='Observed colour · original selections' if frame=='observed' else 'Rest-frame SDSS colour at zref=0.1 · audited selections'
            method_label='Highest-probability class' if method=='hard' else 'Probability-weighted inferred classes'
            fig.suptitle(f'Loa colour–mass diagrams by inferred environment\n{frame_label} · {method_label}',fontsize=15,y=.97)
            fig.text(.065,.035,'Horizontal cuts: 0.75 and 1.5. Shared axes and count scale. Raw distributions; no mass/redshift matching or volume normalization.\nEnvironment predictions are held fixed. Rest-frame colours retain internal attenuation; these are not dust-free or quenching labels.',fontsize=10)
            save(fig,f'colour_mass_{frame}_{method}')
    fig,axs=plt.subplots(2,3,figsize=(17,9),sharex=True,sharey=True)
    fig.subplots_adjust(left=.065,right=.91,bottom=.12,top=.87,wspace=.13,hspace=.23)
    norm=LogNorm(vmin=1,vmax=max(v.max() for k,v in hist.items() if k.startswith('paired')))
    for i,frame in enumerate(['observed','rest']):
        for j,(key,*_) in enumerate(WINDOWS):
            w=windows[key];tag=f'paired_{frame}_{key}';a=data[key,'paired'];col='GR_OBS' if frame=='observed' else 'GR_REST01';cm=masks(a[col])
            note=f"N = {len(a):,}\nRed {cm['red'].sum():,}; blue {cm['blue'].sum():,}\nOutside axes: {audit[tag]['weight_outside']:.0f}"
            mesh=draw(axs[i,j],hist[tag],w['paired_limit'],f"{w['zlo']:.2f} ≤ z < {w['zhi']:.2f} · log M★ ≥ {w['paired_limit']:.2f}",note,norm)
            if j==0:axs[i,j].set_ylabel('Observed g − r [mag]' if frame=='observed' else 'Rest-frame ⁰·¹(g − r) [mag]')
            if i==1:axs[i,j].set_xlabel('log₁₀ stellar mass / M☉ (CIGALE)')
    fig.colorbar(mesh,cax=fig.add_axes([.93,.18,.012,.64]),label='Galaxies per bin (log scale)')
    fig.suptitle('Observed versus rest-frame colour · identical galaxies within each column',fontsize=15,y=.97)
    fig.text(.065,.035,'Matched valid FastSpecFit colours; common mass threshold ≥ both completeness estimates. Fixed numerical cuts 0.75 and 1.5.\nThe frame and filter system change together: observed Legacy Survey g−r versus rest-frame SDSS ⁰·¹(g−r). No internal dust correction.',fontsize=10)
    save(fig,'colour_mass_observed_rest_paired')
    np.savez_compressed(FIG/'diagram_bins.npz',mass_edges=XE,colour_edges=YE,**hist)
    write('diagram_validation.json',audit)


def trends(windows):
    rows=[]
    for key,*_ in WINDOWS:
        w=windows[key]
        for variant in ['rest_primary','paired_observed','paired_rest']:
            r=json.loads((OUT/f'relationships_{variant}_{key}_cut75.json').read_text())
            fig,axs=plt.subplots(2,4,figsize=(18,8.6))
            for ax,(prop,label) in zip(axs.flat,LABELS.items()):
                if prop=='colour':label='Observed g − r [mag]' if variant=='paired_observed' else 'Rest-frame ⁰·¹(g − r) [mag]'
                notes=[]
                for group,col in [('red','#b52c3c'),('blue','#176bad')]:
                    rec=r[group][prop];adj=rec['controlled']['0.1_soft'];notes.append(f"{group}: {rec['rows']:,} raw / {adj.get('rows',0):,} controlled")
                    raw=rec['raw']['soft'];vs=[('Raw',raw,':','s')] if raw is not None else []
                    if adj['status']=='measured':vs.append(('Controlled',adj['methods']['soft'],'-','o'))
                    for method,v,style,marker in vs:
                        ci=np.array(v['ci95']);x=np.arange(4)+(-.035 if group=='red' else .035)
                        ax.plot(x,v['mean'],marker=marker,ls=style,c=col,label=group+' · '+method);ax.vlines(x,ci[0],ci[1],color=col,alpha=.65)
                ax.set_title(label+'\n'+'\n'.join(notes),fontsize=9);ax.set_xticks(range(4),ENV,rotation=20);ax.grid(alpha=.2)
            axs[0,0].legend(fontsize=8)
            limit=w['rest_limit'] if variant=='rest_primary' else w['paired_limit']
            frame='Observed g−r' if variant=='paired_observed' else 'Rest-frame SDSS ⁰·¹(g−r)'
            sample='Rest-colour completeness sample' if variant=='rest_primary' else 'Same-galaxy observed/rest comparison'
            fig.suptitle(f"{frame} · {w['zlo']:.2f} ≤ z < {w['zhi']:.2f} · log₁₀ M★ ≥ {limit:.2f}\n{sample} · red: 0.75 < colour < 1.5; blue: colour < 0.75\nSolid: within-colour mass/z control (mass panel: z only); dotted: raw · conditional 95% sky-block intervals",fontsize=12,y=.98)
            fig.subplots_adjust(left=.05,right=.99,bottom=.105,top=.79,hspace=.65,wspace=.29)
            fig.text(.5,.02,'Inferred environment probabilities; red/blue curves have different support. Colour is not a quenching label. No causal or calibration claim.',ha='center',fontsize=10)
            save(fig,f'properties_{variant}_{key}')
        # Compact sensitivity plot for the requested red/blue environment trends.
        fig,axs=plt.subplots(2,3,figsize=(15,11))
        for i,(group,col) in enumerate([('red','#b52c3c'),('blue','#176bad')]):
            for j,prop in enumerate(['quenched','logssfr','colour']):
                ax=axs[i,j];factor=100 if prop=='quenched' else 1
                specs=[(cut,f'{width}_{support}',f'{cut:.2f} / {width:.2f} / {support}')
                       for cut in [.75,.70,.80]
                       for width,support in [(.1,'soft'),(.05,'soft'),(.1,'joint'),(.05,'joint')]]
                for y,(cut,mode,name) in enumerate(specs):
                    d=json.loads((OUT/f'relationships_rest_primary_{key}_cut{round(cut*100):02d}.json').read_text())[group][prop]['controlled'][mode]
                    if d['status']=='measured':
                        v=d['methods']['soft'];lo,hi=np.array(v['contrast_ci95'])*factor
                        ax.hlines(y,lo,hi,color=col,lw=1.5);ax.plot(v['contrast']*factor,y,'o',color=col)
                    else:ax.text(.02,y,'No common support',transform=ax.get_yaxis_transform(),fontsize=8)
                ax.axvline(0,ls=':',color='gray');ax.set_yticks(range(len(specs)),[s[2] for s in specs],fontsize=8);ax.invert_yaxis()
                label={'quenched':'Δ quenched [percentage points]','logssfr':'Δ log sSFR [dex]','colour':'Δ rest-frame colour [mag]'}[prop]
                ax.set_title(group.capitalize()+' · '+label,fontsize=11);ax.grid(axis='x',alpha=.2)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        fig.suptitle(f"Rest-frame cluster − void contrasts · {w['zlo']:.2f} ≤ z < {w['zhi']:.2f}\nSensitivity: colour boundary / mass-cell width [dex] / support mode · conditional 95% intervals",fontsize=13)
        fig.tight_layout(rect=[0,.05,1,.90]);fig.text(.5,.015,'Joint means soft + hard common support; plotted estimate remains soft. All variants use the same audited mass cut.',ha='center',fontsize=9)
        save(fig,f'contrasts_rest_sensitivity_{key}')
    for path in sorted(OUT.glob('relationships_*.json')):
        r=json.loads(path.read_text())
        for group,props in r.items():
            for prop,rec in props.items():
                for mode,methods in [('raw',rec['raw'])]+[(k,d.get('methods',{})) for k,d in rec['controlled'].items()]:
                    for method,v in methods.items():
                        if v is None:continue
                        rows.append(dict(source=path.name,group=group,property=prop,mode=mode,method=method,
                            void=v['mean'][0],wall=v['mean'][1],filament=v['mean'][2],cluster=v['mean'][3],
                            cluster_minus_void=v['contrast'],ci95_low=v['contrast_ci95'][0],ci95_high=v['contrast_ci95'][1]))
    with (FIG/'curves_and_contrasts.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    global FIG
    ap=argparse.ArgumentParser();ap.add_argument('--output-name',default='presentation');args=ap.parse_args()
    if Path(args.output_name).name!=args.output_name:raise ValueError('Output must be one directory name')
    FIG=OUT/args.output_name
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):raise RuntimeError('Approved compute allocation required')
    FIG.mkdir(exist_ok=False)
    windows=json.loads((OUT/'windows.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    diagrams(windows);trends(windows)
    d=json.loads((OUT/'fixed_mass_colour_fractions.json').read_text())
    fig,ax=plt.subplots(figsize=(9,5))
    for frame,col,marker in [('observed','#555555','s'),('rest','#176bad','o')]:
        rows=[b for b in d['bins'] if b['rows']]
        ax.plot([(b['zlo']+b['zhi'])/2 for b in rows],
            [100*b[frame]['blue_fraction_all'] for b in rows],color=col,marker=marker,label=frame.capitalize())
    ax.set(xlabel='Redshift',ylabel='Blue fraction of matched sample [%]',
        title=f"Fixed numerical colour < 0.75 · identical galaxies in each redshift bin\nFixed log₁₀ M★ ≥ {d['logmass']:.2f} across 0.15 ≤ z < 0.30",ylim=(0,None))
    ax.legend();ax.grid(alpha=.2)
    fig.text(.5,.02,'Descriptive counts; denominator includes colours outside the red/blue cuts. No volume weighting.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.06,1,1]);save(fig,'blue_fraction_fixed_mass')
    write('presentation_manifest.json' if FIG.name=='presentation' else FIG.name+'_manifest.json',dict(job=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),code_sha256=sha(Path(__file__)),
        analysis_sha256=sha(OUT/'analysis_manifest.json'),outputs={p.name:sha(p) for p in FIG.iterdir() if p.is_file()}))
    print('Rendering complete',str(FIG),flush=True)

if __name__=='__main__':main()
