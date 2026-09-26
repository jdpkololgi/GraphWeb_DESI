"""Plot upstream raw-parent upper bound and magnitude sensitivity."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parents[2]/'docs/figures/p12a_selection_20260925'
def main():
 c=json.loads((OUT/'COUNTS.json').read_text());u=json.loads((OUT/'UPSTREAM.json').read_text());z=(np.array(c['edges'][1:])+np.array(c['edges'][:-1]))/2;loa=np.array(c['loa']);mock=np.array(c['mock'])[-1];raw=np.array(u['counts']['raw_rlt19.5'])[:,:40];fig,ax=plt.subplots(2,2,figsize=(12,8),sharex=True,layout='constrained')
 for cap,name in enumerate(['SGC','NGC']):
  a=ax[0,cap]
  for h,ls,col,label in [(loa,'-','black','DESI successes'),(mock,'-','#0072b2','ph006 processed successes'),(raw,'--','#d55e00','ph006 raw r < 19.5 (before assignment)')]:a.plot(z,h[cap],ls=ls,color=col,label=label)
  a.set(yscale='log',ylabel='Galaxies per Δz = 0.01',title=name);a.legend(fontsize=8)
  a=ax[1,cap]
  for limit,col in [(19.5,'#0072b2'),(19.54,'#d55e00'),(19.6,'#009e73'),(19.7,'#cc79a7')]:
   h=np.array(u['counts'][f'raw_rlt{limit}'])[:,:40];a.plot(z,loa[cap]/h[cap],color=col,label=f'DESI / raw r < {limit}')
  a.axhline(1,color='grey',ls=':');a.set(xlabel='Observed redshift',ylabel='DESI successes / raw parent counts');a.legend(fontsize=8)
 fig.suptitle('The high-redshift deficit is already present in the raw bright mock population',fontsize=15);fig.supxlabel('Same sky intersection as the main count plot. Raw curves precede targeting/fibre assignment and are not observed mocks.\nMagnitude changes are sensitivity tests only: apparent-magnitude definitions have not been shown equivalent.',fontsize=9)
 for ext in ['png','pdf']:fig.savefig(OUT/f'05_upstream_counts.{ext}',dpi=170,bbox_inches='tight')
 print('high z DESI / raw',loa[:,-10:].sum(1)/raw[:,-10:].sum(1))
if __name__=='__main__':main()
