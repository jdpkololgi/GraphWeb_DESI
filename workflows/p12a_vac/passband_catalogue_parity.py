"""Verify the FastSpecFit pilot's observed photometry against actual VAC source."""
from pathlib import Path
import json,hashlib
import numpy as np,fitsio
from same_galaxy_passbands import OUT,CAT
import loa_trial as t

def main():
 a=np.load(OUT/'SAME_GALAXY_SAMPLE.npz');cols=['TARGETID','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R']
 with fitsio.FITS(CAT) as f:m=f['METADATA'].read(rows=a['catalogue_rows'],columns=cols)
 ids=m['TARGETID'];p=t.read(Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1/INPUTS_READY.json'))['source']['path'];chunks=[]
 with fitsio.FITS(p) as f:
  for start in range(0,f[1].get_nrows(),250000):
   d=f[1].read(rows=np.arange(start,min(start+250000,f[1].get_nrows())),columns=cols);chunks.append(d[np.isin(d['TARGETID'],ids)])
 d=np.concatenate(chunks);d.sort(order='TARGETID');m.sort(order='TARGETID');assert np.array_equal(m['TARGETID'],d['TARGETID'])
 result=dict(source=p,matched=len(d),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),columns={c:dict(equal=int(np.sum(m[c]==d[c])),max_abs_difference=float(np.max(abs(m[c]-d[c])))) for c in cols[1:]})
 np.savez_compressed(OUT/'PARITY_ROWS.npz',fastspecfit=m,loa=d)
 gr=lambda x:-2.5*np.log10((x['FLUX_G']/x['MW_TRANSMISSION_G'])/(x['FLUX_R']/x['MW_TRANSMISSION_R']))
 result['correct_fsf_minus_loa_dered_gr']=dict(median=float(np.median((-2.5*np.log10(m['FLUX_G']/m['FLUX_R']))-gr(d))),max_abs=float(np.max(abs((-2.5*np.log10(m['FLUX_G']/m['FLUX_R']))-gr(d)))))
 result['correct_flux_relative_residual']={c:float(np.max(abs(m['FLUX_'+c]/(d['FLUX_'+c]/m['MW_TRANSMISSION_'+c])-1))) for c in ['G','R']}
 (OUT/'CATALOGUE_PARITY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
