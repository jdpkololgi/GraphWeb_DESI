"""Plot reference-catalogue counts without conflating BAO and full BRIGHT."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(root):
 d=json.loads((root/'COMPLETE.json').read_text());cat=d['catalogues'];z=np.asarray(d['z_edges']);mid=(z[:-1]+z[1:])/2
 old=np.load(Path(__file__).resolve().parents[2]/'docs/evidence/p12a_joint_flags_20260926/HISTOGRAMS.npz')
 v11=np.asarray(cat['loa_v1p1_full_bright']['count']);v21=old['combinations'][...,0];mock=old['mock_observed'].sum(-1)
 summary={'full_bright':{'v1p1_shell_counts':v11[:,5:45].reshape(2,4,10).sum(-1).tolist(),'v2p1_shell_counts':v21.reshape(2,4,10).sum(-1).tolist(),'v1p1_to_old_mock':(v11[:,5:45].reshape(2,4,10).sum(-1)/mock.reshape(2,4,10).sum(-1)).tolist(),'v1p1_to_v2p1':(v11[:,5:45].reshape(2,4,10).sum(-1)/v21.reshape(2,4,10).sum(-1)).tolist()},'bao':{}}
 fig,ax=plt.subplots(2,3,figsize=(15,8),layout='constrained')
 for i,c in enumerate(['SGC','NGC']):
  ax[i,0].plot(mid[5:45],v21[i]/mock[i],label='Our Loa v2.1 / original mock')
  ax[i,0].plot(mid[5:45],v11[i,5:45]/mock[i],ls='--',label='Loa v1.1 / original mock')
  ax[i,0].set_title(c+' — full BGS-BRIGHT');ax[i,0].set_ylabel('Data / mock counts')
  for j,field in enumerate(['count','weighted'],1):
   obs=np.asarray(cat['loa_v1p1_bao_'+c][field])[i]
   for version,tracer,color in [('kibo-v1','BGS_ANY-02','tab:purple'),('loa-v1','BGS_BRIGHT-02','tab:orange')]:
    vals=[];phases=[]
    for phase in d['phases']:
     item=cat[f'ph{phase:03d}_{version}_{tracer}_{c}']
     if not item.get('missing'):vals.append(item[field][i]);phases.append(phase)
    if not vals:continue
    vals=np.asarray(vals);rat=np.divide(vals,obs,out=np.full_like(vals,np.nan,dtype=float),where=obs>0);mean=rat.mean(0);std=rat.std(0,ddof=1)
    ax[i,j].plot(mid,mean,label=version+' '+tracer,color=color)
    ax[i,j].fill_between(mid,mean-std,mean+std,color=color,alpha=.2)
    den=obs.reshape(5,10).sum(-1);num=vals.reshape(-1,5,10).sum(-1)
    sh=np.divide(num,den,out=np.full_like(num,np.nan,dtype=float),where=den>0)
    item={'phases':phases,'ratios_by_phase_shell':[[float(x) if np.isfinite(x) else None for x in row] for row in sh], 'data_shell_counts_or_weights':den.tolist()}
    for label,sl in [('0p1_0p4',slice(0,30)),('0p4_0p5',slice(30,40))]:
     v=vals[:,sl].sum(-1)/obs[sl].sum();item[label]={'mock_over_data_mean':float(v.mean()),'mock_over_data_sd':float(v.std(ddof=1)),'mock_over_data_min':float(v.min()),'mock_over_data_max':float(v.max())}
    summary['bao'][f'{c}_{version}_{field}']=item
   ax[i,j].set_title(c+' — BAO-selected '+('raw counts' if field=='count' else 'stored WEIGHT sums'))
   ax[i,j].set_ylabel('Mock / Loa v1.1 BAO sample');ax[i,j].axvspan(.4,.6,color='grey',alpha=.13,label='Outside BAO validation range' if i==0 and j==1 else None)
   ax[i,j].set_xlim(.1,.5)
  for j in range(3):ax[i,j].axhline(1,color='black',lw=.8);ax[i,j].set_xlabel('Redshift');ax[i,j].grid(alpha=.2)
 for j in range(3):ax[0,j].legend(fontsize=7)
 fig.suptitle('Loa version change versus BAO mock selection: different scientific samples\nFixed common sky; bands = scatter across exposed ph002–006, not error on the mean',fontsize=13)
 fig.savefig(root/'reference_comparison.png',dpi=160);fig.savefig(root/'reference_comparison.pdf');plt.close(fig)
 (root/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();main(a.root)
