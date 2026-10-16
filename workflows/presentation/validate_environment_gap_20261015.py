"""Independent direct aggregation and source parity checks for new diagnostics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np,json,fitsio
H=Path('/global/u2/d/dkololgi');D=H/'GraphWeb_DESI/docs/figures/environment_gap_three_tests_20261015';R=H/'fmsu_tokyo_20261010';d=np.load(R/'environment_gap_three_tests_20261015/uncertainty_sample.npz');b=json.loads((D/'test3.json').read_text());ex=json.loads((D/'extraction.json').read_text())
checks=[]
def check(name,ok):
    assert ok,name
    checks.append(name)
with fitsio.FITS(ex['source']) as f:
    ii=np.linspace(0,len(d['TARGETID'])-1,257).astype(int)
    a=f[1].read(rows=d['source_row'][ii],columns=['TARGETID','Z','MASS_CG_15','SFR_CG_15','P_WEB','EIGENVALUE_Q16','EIGENVALUE_Q84'])
    for k in a.dtype.names:check('source_'+k,np.array_equal(a[k],d[k][ii],equal_nan=True))
m=d['MASS_CG_15'];sf=d['SFR_CG_15'];mass=np.log10(np.maximum(m,1e-30));z=d['Z'];p=d['P_WEB'];eig=d['EIGENVALUE_MEAN'];width=(d['EIGENVALUE_Q84'][:,0]-d['EIGENVALUE_Q16'][:,0])/2
valid=d['CIGALE_MATCHED'].astype(bool)&~d['CIGALE_AMBIGUOUS'].astype(bool)&(np.char.strip(d['CIGALE_SPECTYPE'].astype(str))=='GALAXY')&np.isfinite(m)&(m>0)&np.isfinite(sf)&(sf>=0)&(mass>=9)&(mass<12)&np.isfinite(eig).all(1)&np.isfinite(width)&(width>0)&np.isfinite(p).all(1)&(p>=0).all(1)&(np.abs(p.sum(1)-1)<1e-4)
cells=np.floor((np.where(np.isfinite(mass),mass,0)-9)/.2+1e-7).astype(int)+100*np.floor((z-.15)/.05+1e-7).astype(int)
low=np.zeros(len(z),bool)
for c in np.unique(cells[valid]):
    ii=np.flatnonzero(valid&(cells==c));order=np.lexsort((d['TARGETID'][ii],width[ii]));low[ii[order[:len(ii)//2]]]=True
q=(sf/np.maximum(m,1e-30)<1e-11).astype(float);maxerr=0
for ri,r in enumerate(b['windows']):
    base=valid&(z>=r['z'][0])&(z<r['z'][1])
    if r['boundary_over_20_mpc']:base&=d['BOUNDARY_MPC']>20
    common=base&np.isin(cells,r['cells'])
    check(f'window_{ri}_N',int(common.sum())==r['common_parent_rows'])
    for arm in r['arms']:
        sel=common&low if arm['low_width'] else common
        vals=np.zeros(4)
        for c in r['cells']:
            cell=sel&(cells==c);w=p[cell];yy=q[cell];vals+=((common&(cells==c)).sum()/common.sum())*np.average(np.tile(yy[:,None],(1,4)),axis=0,weights=w)
        err=np.max(np.abs(vals-arm['class_quenched']));maxerr=max(maxerr,float(err));check(f'window_{ri}_arm_{arm["low_width"]}_means',err<1e-10)
# Ground-truth full-population values must reproduce the earlier independent experiment.
a=json.loads((D/'test1.json').read_text())
for rad,expected in [(2,.4949941019138592),(7,.3272568551096608)]:
    arm=next(r for r in a['arms'] if r['radius']==rad and not r['high_mass'] and r['estimator']=='mass_controlled_trace_quartiles')
    check(f'TNG_R{rad}_prior_parity',abs(arm['contrast'][0]-expected)<1e-12)
check('test1_heldout_N',a['rows']==66691)
check('valid_count',valid.sum()==b['valid_rows'])
check('ordered_lambda',np.all(np.diff(eig[valid],axis=1)>=-1e-6))
v=json.loads((D/'VALIDATION.json').read_text());v.update(independent_checks=checks,independent_check_count=len(checks),maximum_direct_mean_error=maxerr,source_rows_freshly_checked=257)
(D/'VALIDATION.json').write_text(json.dumps(v,indent=2));print(len(checks),'checks; maximum mean error',maxerr)
