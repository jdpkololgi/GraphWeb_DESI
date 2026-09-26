"""Matched-sky number counts and saved ph006 truth/posterior diagnostics."""
import json,sys,hashlib
from pathlib import Path
import numpy as np, fitsio, healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import loa_trial as t
OUT=Path(__file__).resolve().parents[2]/'docs/figures/p12a_selection_20260925'
ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
def save(name,d): (OUT/name).write_text(json.dumps(d,indent=2)+'\n')
def figsave(fig,name):
 for ext in ['png','pdf']:fig.savefig(OUT/(name+'.'+ext),dpi=170,bbox_inches='tight')
 plt.close(fig)
def counts():
 edges=np.linspace(.15,.55,41);zc=(edges[1:]+edges[:-1])/2
 common=np.load(ROOT/'galaxy_angular_support.npy')
 phases=[f'ph{i:03d}' for i in range(2,7)]
 for phase in phases:common &= np.load(t.B/phase/'p3_fields/angular_support_nside256.npz')['support']
 ra,dec=hp.pix2ang(256,np.arange(len(common)),lonlat=True);pc=t.galactic_cap(ra,dec)
 area=np.array([(common&(pc==c)).sum()*hp.nside2pixarea(256) for c in [0,1]])
 def hist(ra,dec,z,cap):
  pix=hp.ang2pix(256,ra,dec,lonlat=True);ok=common[pix]
  return np.array([np.histogram(z[ok&(cap==c)],edges)[0] for c in [0,1]])
 d=np.load(ROOT/'catalogue.npz');loa=hist(d['ra'],d['dec'],d['z'],d['cap']);mock=[];sources=[]
 for phase in phases:
  m=t.read(t.B/phase/'p1_canonical/manifest.json');h=np.zeros((2,40),dtype='i8')
  with fitsio.FITS(m['parent']) as f:
   for start in range(0,f[1].get_nrows(),250000):
    x=f[1][start:min(start+250000,f[1].get_nrows())];k=(x['ZWARN']==0)&np.isfinite(x['Z']);x=x[k]
    h+=hist(x['RA'],x['DEC'],x['Z'],t.galactic_cap(x['RA'],x['DEC']))
  mock.append(h);sources.append(dict(phase=phase,path=m['parent'],sha256=m['parent_sha256']));print(phase,h.sum(axis=1).tolist(),flush=True)
 mock=np.array(mock);selection=t.read(t.read(ROOT/'INPUTS_READY.json')['selection_path']);expected=np.zeros((2,40))
 for c,name in enumerate(['SGC','NGC']):
  curve=selection['rotations']['0']['caps'][name]
  for j,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
   zg=np.linspace(a,b,101);r=np.interp(zg,selection['cosmology']['redshift_grid'],selection['cosmology']['radius_grid_mpc']);n=np.interp(zg,curve['grid_z'],curve['ntilde']);expected[c,j]=area[c]*np.trapz(n*r*r,r)
 mean=mock.mean(axis=0);fig,ax=plt.subplots(2,2,figsize=(12,8),sharex=True,layout='constrained')
 for c,name in enumerate(['SGC','NGC']):
  a=ax[0,c];a.plot(zc,loa[c],color='black',label='DESI Loa');a.plot(zc,mean[c],color='#0072b2',label='Abacus ph002–006 mean');a.fill_between(zc,mock[:,c].min(0),mock[:,c].max(0),color='#0072b2',alpha=.25,label='Five-phase range');a.plot(zc,expected[c],ls='--',color='#d55e00',label='Frozen expected (geometric)');a.set(title=f'{name}: {area[c]*(180/np.pi)**2:,.0f} deg² common sky',ylabel='Galaxies per Δz = 0.01');a.legend(fontsize=9)
  a=ax[1,c];a.plot(zc,loa[c]/mean[c],color='black',label='DESI / mean Abacus');a.fill_between(zc,loa[c]/mock[:,c].max(0),loa[c]/mock[:,c].min(0),color='#0072b2',alpha=.25);a.plot(zc,loa[c]/expected[c],color='#d55e00',ls='--',label='DESI / frozen expected');a.axhline(1,color='grey',ls=':');a.set(xlabel='Observed redshift',ylabel='Count ratio');a.legend(fontsize=9)
 fig.suptitle('The excess grows with redshift; it is not a survey-wide factor of two',fontsize=16)
 fig.supxlabel('Identical NSIDE=256 occupied-sky intersection for DESI and all five mocks. Phase range is not a systematic-error band.\nExpected curve integrates frozen ñ(z) over geometric volume; it does not reproduce voxel apodization.',fontsize=9)
 figsave(fig,'01_matched_counts')
 coarse=lambda x:x.reshape(*x.shape[:-1],4,10).sum(-1)
 report=dict(edges=edges.tolist(),area_deg2=(area*(180/np.pi)**2).tolist(),loa=loa.tolist(),mock=mock.tolist(),expected=expected.tolist(),sources=sources,coarse_loa_over_mock=(coarse(loa)/coarse(mean)).tolist(),coarse_loa_over_expected=(coarse(loa)/coarse(expected)).tolist(),total_loa_over_mock=float(loa.sum()/mean.sum()),script_sha256=t.digest(__file__))
 save('COUNTS.json',report);print('COUNTS',report['coarse_loa_over_mock'],report['total_loa_over_mock'],flush=True)
def phase():
 import torch
 from workflows.sbi.p12_train_base_response_fmpe import theta_to_eigenvalues
 a=t.read(t.C/'posterior/calibration_audit/P12A_CALIBRATION_AUDIT.json');ix=np.load(a['provenance']['evaluation_index']);samples=np.load(a['provenance']['samples'],mmap_mode='r');ck=torch.load(t.C/'posterior/fmpe_estimator.pt',map_location='cpu',weights_only=False)
 d=np.load(t.C/'dataset/ph006_selection_sample.npz');truth=d['truth_eigenvalues'][ix];base=d['base_prediction_eigenvalues'][ix];z=d['context'][ix,3];w=d['natural_weight'][ix];p=[];q=[]
 for start in range(0,len(ix),512):
  eig=theta_to_eigenvalues(samples[start:start+512]*np.asarray(ck['theta_std'])+np.asarray(ck['theta_mean']));cls=(eig>.2).sum(-1);p.append(np.stack([(cls==j).mean(1) for j in range(4)],1));q.append(np.quantile(eig,[.05,.16,.5,.84,.95],axis=1).transpose(1,0,2))
 p=np.concatenate(p);q=np.concatenate(q);tc=(truth>.2).sum(1);bc=(base>.2).sum(1);bins=np.linspace(.15,.55,17);zc=(bins[1:]+bins[:-1])/2;rows=[]
 for lo,hi in zip(bins[:-1],bins[1:]):
  k=(z>=lo)&(z<hi);ww=w[k];avg=lambda v:np.average(v[k],axis=0,weights=ww)
  rows.append(dict(n=int(k.sum()),truth=avg(np.eye(4)[tc]).tolist(),base=avg(np.eye(4)[bc]).tolist(),posterior=avg(p).tolist(),width=avg(q[:,3]-q[:,1]).tolist(),coverage68=avg((truth>=q[:,1])&(truth<=q[:,3])).tolist(),coverage90=avg((truth>=q[:,0])&(truth<=q[:,4])).tolist()))
 desi=t.read(OUT.parent/'p12a_loa_full_20260925/PLOT_SUMMARY.json');cols=['#4477aa','#ccbb44','#228833','#aa3377'];names=['Void','Sheet','Filament','Knot'];fig,ax=plt.subplots(2,2,figsize=(12,8),sharex=True,layout='constrained')
 for j,a in enumerate(ax.flat):
  for key,ls,label in [('truth','-','ph006 truth'),('base',':','ph006 encoder classes'),('posterior','--','ph006 mean posterior')]:a.plot(zc,[r[key][j] for r in rows],ls=ls,color=cols[j],label=label,lw=2)
  a.plot(desi['redshift_centres'],np.array(desi['class_probabilities_by_redshift'])[:,j],color='black',label='DESI mean posterior',lw=1.5);a.set(title=names[j],ylabel='Fraction / mean probability',xlabel='Observed redshift',ylim=(0,.8));a.legend(fontsize=8)
 fig.suptitle('Known mock truth separates population trends from inference effects',fontsize=16);fig.supxlabel('ph006: same 50,000 saved evaluation rows and 512 draws, natural sampling weights; previously exposed development phase.\nDESI: all supported galaxies. Different populations; posterior means are descriptive, not population-PDF estimates.',fontsize=9);figsave(fig,'02_truth_encoder_posterior')
 fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for j in range(3):
  ax[0].plot(zc,[r['coverage68'][j] for r in rows],label=f'λ{j+1} 68%');ax[1].plot(zc,[r['coverage90'][j] for r in rows],label=f'λ{j+1} 90%')
 for a,nom in zip(ax,[.68,.9]):a.axhline(nom,color='black',ls=':');a.set(xlabel='Observed redshift',ylabel='Weighted empirical coverage');a.legend()
 fig.suptitle('ph006 fine-redshift coverage — descriptive development diagnostic');figsave(fig,'03_mock_coverage');save('PH006.json',dict(centres=zc.tolist(),bins=rows,rows=len(ix),audit_path=str(t.C/'posterior/calibration_audit/P12A_CALIBRATION_AUDIT.json'),natural_weighted=True));print('PHASE',rows[0],rows[-1],flush=True)
if __name__=='__main__':
 OUT.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.15})
 {'counts':counts,'phase':phase}[sys.argv[1]]()
