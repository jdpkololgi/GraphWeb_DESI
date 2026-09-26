"""TARGETID join to DR9/1.1.1 and BGS replay including Gaia/SGA."""
from pathlib import Path
import ast,json,logging,os
import numpy as np,fitsio,healpy as hp
import loa_trial as t
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_closure_20260926'
TARGS=Path('/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright')
def main():
 # Use the historical function bodies, not a rewritten selection approximation.
 ns={'np':np,'log':logging.getLogger(__name__)}
 names={'isBGS','notinBGS_mask','isBGS_colors','isBGS_sga','_check_BGS_targtype'}
 source=OUT/'desitarget_1p1p1_cuts.txt';tree=ast.parse(source.read_text());body=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in names];assert len(body)==5
 # Exact BGS geometry mask in tag1.1.1 is BRIGHT(bit1), CLUSTER(bit13).
 ns['imaging_mask']=lambda maskbits,bgsmask=False: (maskbits&((1<<1)|(1<<13)))==0
 exec(compile(ast.Module(body=body,type_ignores=[]),str(source),'exec'),ns)
 root=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1');cat=np.load(root/'catalogue.npz');take=cat['shell']>=0;ids=cat['targetid'][take];z=cat['z'][take];cap=cat['cap'][take];ra=cat['ra'][take];dec=cat['dec'][take];order=np.argsort(ids);ids=ids[order];z=z[order];cap=cap[order];ra=ra[order];dec=dec[order];shell=np.searchsorted([.15,.25,.35,.45,.55],z,side='right')-1;b=cap*4+shell
 hpix=np.unique(hp.ang2pix(8,ra,dec,nest=True,lonlat=True));seen=np.zeros(len(ids),bool);passed=np.zeros(len(ids),bool);sga=np.zeros(len(ids),bool);stored=np.zeros(len(ids),bool);gaiapass=np.zeros(len(ids),bool);phot=np.zeros((len(ids),3),dtype='f4');run=[];allbit=0;replaydiff=0
 for n,pix in enumerate(hpix):
  path=TARGS/f'targets-bright-hp-{pix}.fits'
  if not path.exists():continue
  with fitsio.FITS(path) as f:d=f[1][:];h=f[1].read_header()
  pos=np.searchsorted(ids,d['TARGETID']);valid=pos<len(ids);valid[valid]&=ids[pos[valid]]==d['TARGETID'][valid];d=d[valid];pos=pos[valid]
  if not len(d):continue
  assert not seen[pos].any();seen[pos]=True
  flux={c:d['FLUX_'+c]/d['MW_TRANSMISSION_'+c] for c in ['G','R','Z','W1']};gaia=d['GAIA_PHOT_G_MEAN_MAG'];grr=gaia-22.5+2.5*np.log10(np.maximum(d['FLUX_R'],1e-16));gaiapass[pos]=(grr>.6)|(gaia==0);sga[pos]=np.char.startswith(d['REF_CAT'].astype('U'),'L');stored[pos]=(d['BGS_TARGET']&2)!=0
  kw=dict(rfiberflux=d['FIBERFLUX_R']/d['MW_TRANSMISSION_R'],rfibertotflux=d['FIBERTOTFLUX_R'],gflux=flux['G'],rflux=flux['R'],zflux=flux['Z'],w1flux=flux['W1'],gnobs=d['NOBS_G'],rnobs=d['NOBS_R'],znobs=d['NOBS_Z'],gfluxivar=d['FLUX_IVAR_G'],rfluxivar=d['FLUX_IVAR_R'],zfluxivar=d['FLUX_IVAR_Z'],maskbits=d['MASKBITS'],Grr=grr,refcat=d['REF_CAT'],gaiagmag=gaia,targtype='bright')
  north=np.char.strip(d['PHOTSYS'].astype('U'))=='N';got=np.where(north,ns['isBGS'](**kw,south=False),ns['isBGS'](**kw,south=True));passed[pos]=got
  for j,c in enumerate(['G','R','Z']):phot[pos,j]=d['FLUX_'+c]
  run.append(dict(path=str(path),rows=len(d),desitarget=h.get('DEPVER10'),photcat=h.get('DEPVER13'),command=h.get('CMDLINE')))
  if n%40==0:print('targetfiles',n,len(hpix),'matched',int(seen.sum()),flush=True)
 # Check original versus LSS fluxes by exact TARGETID. Output differences only.
 mismatch=np.zeros(3,dtype='i8');maxdiff=np.zeros(3);sourcepath=t.read(root/'INPUTS_READY.json')['source']['path']
 with fitsio.FITS(sourcepath) as f:
  for start in range(0,f[1].get_nrows(),250000):
   d=f[1][start:min(start+250000,f[1].get_nrows())];pos=np.searchsorted(ids,d['TARGETID']);k=pos<len(ids);k[k]&=ids[pos[k]]==d['TARGETID'][k];d=d[k];pos=pos[k];k=seen[pos];pos=pos[k];d=d[k]
   for j,c in enumerate(['G','R','Z']):dif=abs(phot[pos,j]-d['FLUX_'+c]);mismatch[j]+=(dif!=0).sum();maxdiff[j]=max(maxdiff[j],float(dif.max(initial=0)))
 out=dict(job=os.environ.get('SLURM_JOB_ID'),rows=len(ids),matched=int(seen.sum()),unmatched=int((~seen).sum()),matched_pass=int((seen&passed).sum()),matched_fail=int((seen&~passed).sum()),stored_bright_mismatch=int((seen&~stored).sum()),replay_stored_mismatch=int((seen&(passed!=stored)).sum()),gaia_fail=int((seen&~gaiapass).sum()),sga_rows=int((seen&sga).sum()),gaia_fail_recovered_sga=int((seen&~gaiapass&sga&passed).sum()),fail_by_shell=np.bincount(b[seen&~passed],minlength=8).tolist(),flux_changed_rows=mismatch.tolist(),flux_max_abs_diff=maxdiff.tolist(),files=run,script_sha256=t.digest(__file__),historical_cuts_sha256=t.digest(source),historical_geomask_sha256=t.digest(OUT/'desitarget_1p1p1_geomask.txt'),note='Join and replay on selected Loa science galaxies. Does not reconstruct unselected galaxy completeness or prove mock equivalence. Data release header and sampled FA header pin DR9/1.1.1. Remaining unmatched IDs must not be presumed to share this provenance.')
 t.save(OUT/'TARGET_REPLAY.json',out);print({k:v for k,v in out.items() if k!='files'},flush=True)
if __name__=='__main__':main()
