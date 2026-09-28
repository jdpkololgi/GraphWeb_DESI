"""Loa observed and assignment-corrected counts on the existing common sky."""
import json,os,hashlib
from pathlib import Path
import numpy as np, fitsio, healpy as hp
from mock_selection_test import cap
from alignment_parent_screen import MASK
R=Path(__file__).resolve().parents[2]
def main():
 path=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
 mask=np.load(MASK);edges=np.linspace(.15,.55,41);h={};bad=0;selected=0
 cols=['RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE','FLUX_R','MW_TRANSMISSION_R','FRACZ_TILELOCID','FRAC_TLOBS_TILES','WEIGHT_ZFAIL']
 with fitsio.FITS(path) as f:
  for start in range(0,f[1].get_nrows(),250000):
   d=f[1].read(rows=np.arange(start,min(start+250000,f[1].get_nrows())),columns=cols)
   z=d['Z_not4clus'];ok=(z>=.15)&(z<.55)&(d['ZWARN']==0)&(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
   d=d[ok];d=d[mask[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]]
   z=d['Z_not4clus'];c=cap(d['RA'],d['DEC']);a=d['FRACZ_TILELOCID'];b=d['FRAC_TLOBS_TILES'];valid=np.isfinite(a*b)&(a>0)&(b>0)&(a<=1)&(b<=1)
   bad+=int((~valid).sum());selected+=len(d)
   assert np.all(valid),'invalid completeness among selected rows'
   r=22.5-2.5*np.log10(d['FLUX_R']/d['MW_TRANSMISSION_R'])
   weights={'observed':np.ones(len(d)),'assignment':1/(a*b),'assignment_zfail':d['WEIGHT_ZFAIL']/(a*b)}
   for cut,k in [('targeted',np.ones(len(d),bool)),('uniform_r19p5',(r>=12)&(r<19.5))]:
    for label,w in weights.items():
     assert np.all(np.isfinite(w)&(w>0))
     key=cut+'_'+label;h.setdefault(key,np.zeros((2,40)))
     for i in [0,1]:h[key][i]+=np.histogram(z[k&(c==i)],edges,weights=w[k&(c==i)])[0]
   print(start,flush=True)
 out={'source':str(path),'job':os.environ.get('SLURM_JOB_ID'),'selected':selected,'invalid_completeness':bad,'definition':'1/(FRACZ_TILELOCID*FRAC_TLOBS_TILES), WEIGHT_ZFAIL separate; no FKP or imaging weight; assignment-corrected spectroscopic sample, not imaging-complete truth','shell_counts':{k:v.reshape(2,4,10).sum(-1).tolist() for k,v in h.items()},'fine_counts':{k:v.tolist() for k,v in h.items()},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 for v in ['v0.1','v1']:
  counts=np.array(json.load(open(R/f'docs/evidence/p12a_alignment_execution_20260926/parent/{v}.json'))['counts'])
  out[v+'_ratio']={k:(counts/np.array(n)).tolist() for k,n in out['shell_counts'].items()}
 p=R/'docs/evidence/p12a_loa_completeness_20260928';p.mkdir(exist_ok=True);(p/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
