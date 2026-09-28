"""Compare conditional colour residuals and common-mapping luminosity shapes."""
import sys,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from shared.plot_style import apply_style,finalize_axes
apply_style();palette=plt.rcParams['axes.prop_cycle'].by_key()['color'];colors={'v0.1':palette[0],'v1':palette[1],'loa':palette[2]};root=R/'docs/evidence/p12a_joint_colour_luminosity_20260928'
x=json.loads((root/'RESULTS.json').read_text());h=np.load(root/'HISTOGRAMS.npz')
f,ax=plt.subplots(2,2,figsize=(14,9),layout='constrained')
for i,p in enumerate(['S','N']):
 for v in ['v0.1','v1']:
  rows=[r for r in x['shared_rows'] if r['photsys']==p and r['mock']==v];z=[(r['zlo']+r['zhi'])/2 for r in rows]
  ax[0,i].plot(z,[r['mean_colour_mock_minus_loa'] for r in rows],marker='o',label=v,color=colors[v])
 ax[0,i].set_ylim(-.13,.04)
 ax[0,i].axhline(0,color='white',ls=':',lw=1)
 finalize_axes(ax[0,i],f'PHOTSYS {p}: matched z/r cells','Redshift','Mean Δ(g-r), mock − Loa [mag]')
 m=(h['M_edges'][:-1]+h['M_edges'][1:])/2
 for v in ['loa','v0.1','v1']:
  y=h[v+'_M'][i,30:40].sum(0);ax[1,i].plot(m,y/y.sum(),label=v,color=colors[v])
 ax[1,i].set_xlim(-25.5,-20)
 finalize_axes(ax[1,i],f'PHOTSYS {p}: 0.45 ≤ z < 0.55','Common-mapping M_r - 5 log h','Fraction per 0.1 mag')
 ax[0,i].legend();ax[1,i].legend()
f.suptitle('Exploratory stride-20 sample; raw parents vs completeness-corrected Loa\nColour conditioned on z/r; luminosity uses the same v1 K/E hypothesis, not an independent LF',fontsize=11)
out=R/'docs/figures/p12a_joint_colour_luminosity_20260928';out.mkdir(exist_ok=True)
for ext in ['png','pdf']:f.savefig(out/('joint.'+ext),dpi=160,bbox_inches='tight')
