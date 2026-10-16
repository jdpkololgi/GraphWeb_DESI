"""Final presentation and tabular export from immutable colour analysis results."""
import os
os.environ['MPLBACKEND']='Agg'
import csv,json,socket
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from colour_split import LABELS
from colour_followup import OUT
from run import sha,write_json

def main():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):raise RuntimeError('Compute allocation required')
    out=OUT/'presentation';out.mkdir(exist_ok=False)
    data=json.loads((OUT/'relationships.json').read_text());windows=json.loads((OUT/'windows.json').read_text())
    rows=[]
    for key,groups in data.items():
        fig,axs=plt.subplots(2,4,figsize=(18,8.6),layout=None)
        for ax,(prop,label) in zip(axs.flat,LABELS.items()):
            notes=[]
            for group,col in [('red','#b52c3c'),('blue','#176bad')]:
                rec=groups[group][prop];adj=rec['controlled']['0.1_soft'];notes.append(f"{group}: {rec['rows']:,} raw / {adj.get('rows',0):,} controlled")
                for method,v,style,marker in [('Raw',rec['raw']['soft'],':','s')]+([('Controlled',adj['methods']['soft'],'-','o')] if adj['status']=='measured' else []):
                    ci=np.array(v['ci95']);x=np.arange(4)+(-.035 if group=='red' else .035)
                    ax.plot(x,v['mean'],marker=marker,ls=style,c=col,label=group+' · '+method);ax.vlines(x,ci[0],ci[1],color=col,alpha=.65)
                for mode,record in [('raw',rec['raw'])]+[(k,v.get('methods',{})) for k,v in rec['controlled'].items()]:
                    for method,v in record.items():
                        rows.append({'window':key,'colour':group,'property':prop,'mode':mode,'method':method,'void':v['mean'][0],'wall':v['mean'][1],'filament':v['mean'][2],'cluster':v['mean'][3],'cluster_minus_void':v['contrast'],'ci95_low':v['contrast_ci95'][0],'ci95_high':v['contrast_ci95'][1]})
            ax.set_title(label+'\n'+'\n'.join(notes),fontsize=9);ax.set_xticks(range(4),['Void','Wall','Filament','Cluster'],rotation=20);ax.grid(alpha=.2)
        axs[0,0].legend(fontsize=8);w=windows[key]
        fig.suptitle(f"Loa 0.20 ≤ z < {w['zhi']:.2f} · log₁₀ M★ ≥ {w['logmass']:.2f}\nRed: 0.75 < g−r < 1.5; blue: g−r < 0.75 · inferred environment probabilities\nSolid: mass/z control within each colour (mass panel: z only); dotted: raw · 95% sky-block intervals",fontsize=12,y=.98)
        fig.subplots_adjust(left=.05,right=.99,bottom=.105,top=.79,hspace=.65,wspace=.29)
        fig.text(.5,.02,'Red/blue controlled curves have different support. No tested window passed the predeclared blue-count gate.',ha='center',fontsize=10)
        for suffix in ['png','pdf']:fig.savefig(out/f'properties_{key}.{suffix}',dpi=180)
        plt.close(fig)
    key='z020022';fig,axs=plt.subplots(2,4,figsize=(15,7))
    for ax,(prop,label) in zip(axs.flat,LABELS.items()):
        factor=100 if prop=='quenched' else 1
        for j,(group,col) in enumerate([('red','#b52c3c'),('blue','#176bad')]):
            v=data[key][group][prop]['controlled']['0.1_soft']['methods']['soft'];ci=np.array(v['contrast_ci95'])*factor
            ax.plot(j,v['contrast']*factor,'o',c=col);ax.vlines(j,ci[0],ci[1],color=col,lw=2)
        ax.axhline(0,c='gray',lw=1,ls=':');ax.set_xticks([0,1],['Red','Blue']);ax.set_xlim(-.6,1.6)
        ax.set_title(label.replace('Quenched fraction (sSFR < 10⁻¹¹ yr⁻¹)','Quenched fraction [percentage points]'),fontsize=9);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Inferred cluster − void contrasts · Loa 0.20 ≤ z < 0.22 · log₁₀ M★ ≥ 10.69\nWithin-colour mass/redshift control; mass panel controls z only · conditional 95% sky-block intervals',fontsize=12)
    fig.subplots_adjust(left=.065,right=.99,bottom=.12,top=.8,hspace=.5,wspace=.5)
    fig.text(.5,.025,'Different red/blue support: these are within-population associations, not a matched comparison of colour effects.',ha='center',fontsize=10)
    for suffix in ['png','pdf']:fig.savefig(out/f'contrasts_z020022.{suffix}',dpi=180)
    plt.close(fig)
    with (out/'curves_and_contrasts.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_json(out/'MANIFEST.json',{'analysis_manifest_sha256':sha(OUT/'MANIFEST.json'),'script_sha256':sha(Path(__file__)),'job':os.environ['SLURM_JOB_ID'],'outputs':{p.name:sha(p) for p in out.iterdir() if p.is_file()}})
if __name__=='__main__':main()
