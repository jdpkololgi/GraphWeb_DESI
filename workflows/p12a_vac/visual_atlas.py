"""Read-only spatial atlas and integrity audit of the frozen provisional Loa VAC."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import fitsio
import numpy as np
from astropy.cosmology import Planck18
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

NAMES = ['Void', 'Sheet', 'Filament', 'Knot']
COLORS = ['#62b7ef', '#edca72', '#e97a65', '#ba8ef0']
BG = '#101c2b'
FOOT = ('Provisional environmental inferences · real-DESI coverage unverified\n'
        'Training-mock mismatch grows above z≈0.35, severe at 0.45<z<0.55. Trends are not evidence of physical evolution.\n'
        'Lower-redshift results are not automatically certified. Morphology is not an independent calibration test.')

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def main(receipt_path, out):
    receipt = json.loads(receipt_path.read_text())
    path = Path(receipt['catalogue'])
    before = digest(path)
    assert before == receipt['catalogue_sha256'], 'Original catalogue hash mismatch'
    d = fitsio.read(path)
    s = d['SUPPORTED']; p = d['P_WEB'][s]
    qcols = ['EIGENVALUE_Q05', 'EIGENVALUE_Q16', 'EIGENVALUE_Q50', 'EIGENVALUE_Q84', 'EIGENVALUE_Q95']
    checks = {
        'historical_catalogue_hash_matches': True,
        'row_count_matches': len(d) == receipt['rows'],
        'supported_count_matches': int(s.sum()) == receipt['supported_rows'],
        'unique_targetids': len(np.unique(d['TARGETID'])) == len(d),
        'coordinates_finite_in_range': bool(np.all(np.isfinite(d['RA']) & np.isfinite(d['DEC']) & np.isfinite(d['Z']) & (d['RA']>=0) & (d['RA']<360) & (np.abs(d['DEC'])<=90) & (d['Z']>=.15) & (d['Z']<.55))),
        'support_flag_consistent': bool(np.array_equal(~s, (d['QUALITY'] & 1) != 0)),
        'all_rows_provisional': bool(np.all((d['QUALITY'] & 32) != 0)),
        'probabilities_finite_bounded_normalized': bool(np.isfinite(p).all() and np.all((p>=0)&(p<=1)) and np.allclose(p.sum(axis=1), 1, atol=1e-6)),
        'unsupported_probabilities_null': bool(np.isnan(d['P_WEB'][~s]).all()),
        'quantiles_finite_ordered': bool(all(np.isfinite(d[k][s]).all() and np.all(np.diff(d[k][s],axis=1)>=-1e-6) for k in qcols)),
        'quantiles_monotonic': bool(all(np.all(d[a][s]<=d[b][s]+1e-6) for a,b in zip(qcols[:-1],qcols[1:]))),
        'quality_counts_match': all(int(np.count_nonzero(d['QUALITY']&int(k)))==v for k,v in receipt['quality_bit_counts'].items()),
    }
    assert all(checks.values()), checks
    out.mkdir(parents=True, exist_ok=False)
    # Display geometry chosen in advance, independent of inferred environments.
    zz = np.linspace(.15,.55,10001)
    radius = np.interp(d['Z'], zz, Planck18.comoving_distance(zz).value)
    ang = np.deg2rad(d['RA']-180); dec = np.deg2rad(d['DEC'])
    x = radius*np.cos(dec)*np.sin(ang)
    y = radius*np.cos(dec)*np.cos(ang)
    height = radius*np.sin(dec)
    slab = (d['RA']>=130)&(d['RA']<=230)&(np.abs(height)<10)
    zoom = slab & (np.abs(x)<180)&(y>=650)&(y<=950)
    cls = np.full(len(d), -1); cls[s] = p.argmax(axis=1)
    confidence = np.full(len(d), np.nan); confidence[s] = p.max(axis=1)
    entropy = np.full(len(d), np.nan)
    entropy[s] = -np.sum(p*np.log(np.clip(p,1e-30,1)),axis=1)/np.log(4)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'figure.facecolor':BG,'axes.facecolor':BG,'savefig.facecolor':BG,'text.color':'#edf2f7','axes.labelcolor':'#dce4ef','xtick.color':'#bac8d8','ytick.color':'#bac8d8','axes.edgecolor':'#708399','grid.color':'#53667a','grid.alpha':.2,'axes.spines.top':False,'axes.spines.right':False})
    pdf = PdfPages(out/'LOA_environment_atlas.pdf'); figures=[]
    cmap = ListedColormap(COLORS); norm = BoundaryNorm(np.arange(-.5,4.5),4)
    def save(fig, name):
        fig.text(.5,.015,FOOT,ha='center',va='bottom',fontsize=9,color='#c2cddc')
        fig.savefig(out/(name+'.png'),dpi=190,facecolor=BG)
        pdf.savefig(fig,facecolor=BG); plt.close(fig); figures.append(name+'.png')
    def axes_style(a, mask):
        a.set_aspect('equal'); a.set_xlabel('Transverse coordinate [comoving Mpc]'); a.grid(True,linewidth=.5)
        a.set_xlim(x[mask].min()-4,x[mask].max()+4); a.set_ylim(y[mask].min()-4,y[mask].max()+4)
    def points(a, mask, what, size=2):
        if what=='positions':
            a.scatter(x[mask],y[mask],s=size,c='#cbd5e1',alpha=.65,lw=0,rasterized=True)
        else:
            v=mask&s
            kw = dict(cmap=cmap,norm=norm) if what=='class' else dict(cmap='cividis',vmin=.25,vmax=1)
            im=a.scatter(x[v],y[v],s=size,c=cls[v] if what=='class' else confidence[v],lw=0,rasterized=True,**kw)
            u=mask&~s
            a.scatter(x[u],y[u],s=size,c='#64748b',marker='x',linewidths=.4,rasterized=True)
            return im
    fig,axs=plt.subplots(1,3,figsize=(17,7.5)); fig.subplots_adjust(left=.06,right=.98,top=.84,bottom=.18,wspace=.23)
    fig.suptitle('DESI Loa | a slice through the observed cosmic web',fontsize=23,y=.97)
    fig.text(.5,.875,'Same galaxies, three views · 20 Mpc thick equatorial slab · RA 130°–230° · 0.15 ≤ z < 0.55\nPositions use observed redshifts and Planck18 distances; redshift-space distortions remain.',ha='center',fontsize=12)
    for a,w,title in zip(axs,['positions','class','confidence'],['Observed galaxy positions','Most probable environment','Largest class probability']):
        im=points(a,slab,w,size=1.6); axes_style(a,slab); a.set_title(title,pad=12)
        for z,ls in [(.35,'--'),(.45,':')]:
            r=Planck18.comoving_distance(z).value; t=np.linspace(-np.deg2rad(50),np.deg2rad(50),300)
            a.plot(r*np.sin(t),r*np.cos(t),ls,color='#b8c5d6',lw=.7)
            a.text(0,r+10,f'z = {z}',ha='center',fontsize=8,color='#c2cddc')
    axs[0].set_ylabel('Radial-plane coordinate [comoving Mpc]')
    cax=fig.add_axes([.73,.12,.22,.016]); fig.colorbar(im,cax=cax,orientation='horizontal',label='Posterior concentration, not verified accuracy')
    handles=[Line2D([],[],marker='o',ls='',color=c,label=n) for c,n in zip(COLORS,NAMES)]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.08,.10),ncol=4,frameon=False)
    save(fig,'01_cosmic_web_slice')
    fig,axs=plt.subplots(1,3,figsize=(17,7.5)); fig.subplots_adjust(left=.06,right=.98,top=.78,bottom=.24,wspace=.2)
    fig.suptitle('A closer view | positions, environments and ambiguity',fontsize=22,y=.97)
    fig.text(.5,.895,'Fixed geometric window: transverse ±180 Mpc, radial-plane 650–950 Mpc; same 20 Mpc slab\nNo confidence cut, class balancing, smoothing or invented filament connections.',ha='center',fontsize=11)
    for a,w,title in zip(axs,['positions','class','confidence'],['Observed galaxy positions','Most probable environment','Largest class probability']):
        im=points(a,zoom,w,size=6); axes_style(a,zoom); a.set_title(title)
    axs[0].set_ylabel('Radial-plane coordinate [comoving Mpc]')
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.06,.13),ncol=4,frameon=False)
    cax=fig.add_axes([.73,.16,.22,.018]); fig.colorbar(im,cax=cax,orientation='horizontal',label='Largest class probability')
    save(fig,'02_cosmic_web_closeup')
    fig,axs=plt.subplots(2,2,figsize=(12,12)); fig.subplots_adjust(left=.08,right=.9,top=.86,bottom=.18,hspace=.27,wspace=.22)
    fig.suptitle('The web is inferred probabilistically',fontsize=23,y=.965)
    fig.text(.5,.92,'Same close-up, all four probabilities · bright points have greater probability of that environment\nVoid denotes a galaxy’s local environment, not a map of empty volume.',ha='center',fontsize=11)
    v=zoom&s
    for j,a in enumerate(axs.flat):
        im=a.scatter(x[v],y[v],c=d['P_WEB'][v,j],s=7,cmap='magma',vmin=0,vmax=1,lw=0,rasterized=True)
        axes_style(a,zoom); a.set_title(NAMES[j]); a.set_ylabel('Radial-plane coordinate [Mpc]')
    cax=fig.add_axes([.93,.25,.014,.52]); fig.colorbar(im,cax=cax,label='Per-galaxy posterior probability')
    save(fig,'03_four_environment_probabilities'); pdf.close()
    after=digest(path); checks['catalogue_unchanged_after_plotting']=before==after
    assert checks['catalogue_unchanged_after_plotting']
    report={'catalogue':str(path),'catalogue_sha256':after,'rows':len(d),'supported':int(s.sum()),'unsupported':int((~s).sum()),'checks':checks,'quality_bits':receipt['quality_bits'],'quality_bit_counts':receipt['quality_bit_counts'],'columns':[{ 'name':n,'dtype':str(d.dtype[n])} for n in d.dtype.names],'geometry':{'RA_degrees':[130,230],'absolute_equatorial_height_Mpc_lt':10,'zoom_x_Mpc':[-180,180],'zoom_y_Mpc':[650,950],'cosmology':'Planck18','distance_unit':'Mpc','redshift_space':True,'posterior_dependent_selection':False,'thinning':False},'views':{},'figures':figures,'script_sha256':digest(__file__),'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'limitations':FOOT.split('\n'),'science_release_ready':False}
    for name,m in [('slab',slab),('closeup',zoom)]:
        report['views'][name]={'rows':int(m.sum()),'supported':int((m&s).sum()),'unsupported':int((m&~s).sum()),'z_min':float(d['Z'][m].min()),'z_max':float(d['Z'][m].max()),'argmax_class_counts':[int(np.sum(m&(cls==j))) for j in range(4)],'fraction_max_probability_below_0p5':float(np.mean(confidence[m&s]<.5)),'normalized_entropy_median':float(np.median(entropy[m&s]))}
    (out/'INTEGRITY_AND_VIEWS.json').write_text(json.dumps(report,indent=2)+'\n')
    manifest={p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out/'ARTIFACT_SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'checks':checks,'views':report['views']},indent=2),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('receipt',type=Path); a.add_argument('output',type=Path)
    args=a.parse_args(); main(args.receipt,args.output)
