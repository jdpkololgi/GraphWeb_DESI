"""Frozen-probability off-slab tracer geometry. Not tidal truth or causal inference."""
import argparse,json,os,sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from prepare_parent import digest,rotation,CENTER,RADIUS


def shape(points):
    """Descending covariance spectrum; L,P,S sum to one, not class probabilities."""
    if len(points)<3:return np.full(3,np.nan),np.full(3,np.nan)
    e,v=np.linalg.eigh(np.cov(np.asarray(points).T));e=np.maximum(e,0)
    if e[-1]<=0:return np.full(3,np.nan),np.full(3,np.nan)
    return np.array([(e[2]-e[1])/e[2],(e[1]-e[0])/e[2],e[0]/e[2]]),v[:,0]


def contrast2(points):
    if len(points)<3:return np.nan
    e=np.linalg.eigvalsh(np.cov(points.T))
    return float((e[-1]-e[0])/e.sum()) if e.sum()>0 else np.nan


def blocks_ci(values,positions,seed=17017):
    ok=np.isfinite(values);values=np.asarray(values)[ok];positions=positions[ok]
    if not len(values):return None
    _,b=np.unique(np.floor(positions/30).astype(int),axis=0,return_inverse=True)
    sums=np.bincount(b,weights=values);n=np.bincount(b);rng=np.random.default_rng(seed)
    ix=rng.integers(len(n),size=(512,len(n)));means=sums[ix].sum(1)/n[ix].sum(1)
    return {'mean':float(values.mean()),'ci95':np.quantile(means,[.025,.975]).tolist(),'blocks':len(n),'rows':len(values)}


def summarize(a,mask,position):
    ok=mask&np.isfinite(a).all(1)
    if not ok.any():return {'n':0}
    return {'n':int(ok.sum()),'mean_LPS':a[ok].mean(0).tolist(),'median_LPS':np.median(a[ok],0).tolist(),
       'dominant_counts':np.bincount(a[ok].argmax(1),minlength=3).tolist(),
       'planarity_minus_linearity':blocks_ci((a[:,1]-a[:,0])[ok],position[ok]),
       'sheet_dominant':blocks_ci((a[ok].argmax(1)==1).astype(float),position[ok])}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Compute allocation required')
    if args.out.exists():raise FileExistsError('Output must be new')
    ready=json.loads((args.parent/'READY.json').read_text());p=args.parent/'parent.npz'
    assert digest(p)==ready['parent_sha256'];d=dict(np.load(p));audit=dict(np.load(args.audit))
    xyz=d['xyz_Mpc'];ti=d['original_target_parent_index'];pos=xyz[ti];n=len(ti)
    np.testing.assert_array_equal(d['TARGETID'][ti],audit['TARGETID']);np.testing.assert_array_equal(d['P_WEB'][ti],audit['P_WEB'])
    prob=d['P_WEB'][ti];line=audit['projected_line_like'];low=prob[:,2]<.2
    cohorts={'low_fil_wall':line&low&(prob.argmax(1)==1),'low_fil_void':line&low&(prob.argmax(1)==0),
             'line_fil_argmax':line&(prob.argmax(1)==2),'all_original':np.ones(n,bool),
             'low_fil_wall_pwall_ge05':line&low&(prob[:,1]>=.5)}
    assert cohorts['low_fil_wall'].sum()==297 and cohorts['low_fil_void'].sum()==185
    # Guard catalogue radial truncation and padding. Angular boundary is inherited,
    # not a new exact survey-mask certificate.
    from astropy.cosmology import Planck18
    limits=Planck18.comoving_distance([.15,.55]).value
    radial=np.linalg.norm(pos,axis=1)
    arrays={'TARGETID':d['TARGETID'][ti],'P_WEB':prob,'xyz_Mpc':pos,**{'cohort_'+k:v for k,v in cohorts.items()}}
    rows=[];tree=cKDTree(xyz);radii=ready['config']['radii_Mpc'];rng=np.random.default_rng(17017)
    for ir,r in enumerate(radii):
        eligible=(np.linalg.norm(pos-CENTER,axis=1)+r<RADIUS)&(radial-r>limits[0])&(radial+r<limits[1])&(d['BOUNDARY_MPC'][ti]>r)
        full=np.full((n,3),np.nan);thin=full.copy();matched=full.copy();normal=full.copy();nf=np.zeros(n,int);ns=nf.copy()
        for i,q in enumerate(tree.query_ball_point(pos,r,workers=8)):
            q=np.asarray(q,dtype=int);q=q[q!=ti[i]];s=q[np.abs(xyz[q,2])<10];nf[i]=len(q);ns[i]=len(s)
            if not eligible[i]:continue
            full[i],normal[i]=shape(xyz[q]);thin[i],_=shape(xyz[s])
            if len(s)>=6:
                matched[i]=np.mean([shape(xyz[rng.choice(q,len(s),replace=False)])[0] for _ in range(32)],axis=0)
        for k,v in [('full',full),('slab',thin),('matched',matched),('normal',normal),('nfull',nf),('nslab',ns),('eligible',eligible)]:arrays[f'{k}_{ir}']=v
        for name,cohort in cohorts.items():
            for minimum in [6,8,12]:
                pair=cohort&eligible&(nf>=minimum)&(ns>=minimum)
                row={'radius_Mpc':r,'cohort':name,'min_neighbors':minimum,'cohort_n':int(cohort.sum()),'geometry_eligible':int((cohort&eligible).sum()),
                 'full':summarize(full,cohort&eligible&(nf>=minimum),pos),
                 'paired_full':summarize(full,pair,pos),'paired_slab':summarize(thin,pair,pos),'paired_count_matched':summarize(matched,pair,pos),
                 'paired_delta_PminusL':blocks_ci(((full[:,1]-full[:,0])-(thin[:,1]-thin[:,0]))[pair],pos[pair]),
                 'matched_delta_PminusL':blocks_ci(((matched[:,1]-matched[:,0])-(thin[:,1]-thin[:,0]))[pair],pos[pair]),
                 'median_full_neighbors':float(np.median(nf[cohort])),'median_offslab_fraction':float(np.median(1-ns[cohort]/np.maximum(nf[cohort],1)))}
                rows.append(row)
        print('radius complete',r,flush=True)
    # Independently recover original statistic on exactly the saved tracers.
    oldtree=cKDTree(pos[:,:2]);oldstat=np.full(n,np.nan)
    for i,q in enumerate(oldtree.query_ball_point(pos[:,:2],radii[1])):
        q=np.asarray(q);q=q[q!=i]
        if len(q)>=8:oldstat[i]=contrast2(pos[q,:2])
    valid=audit['projected_shape_valid']
    np.testing.assert_allclose(oldstat[valid],audit['projected_anisotropy'][valid],atol=1e-10,rtol=0)
    arrays['original_projected_shape']=oldstat
    # New parent views intentionally lift original redshift/RA cuts; report transitions.
    # Re-select independent views from parent; frozen shape threshold retained.
    views=[];cut=ready['config']['shape_cut_frozen'];r=radii[1]
    for angle in [0,30,60,90]:
        coords=(xyz-CENTER)@rotation(angle).T;targets=coords[ti]
        for offset in [-5,0,5]:
            for depth in [5,10,20,40,80]:
                select=(abs(coords[:,0])<180)&(abs(coords[:,1])<=225)&(abs(coords[:,2]-offset)<depth/2)
                ix=np.flatnonzero(select);tr=cKDTree(coords[ix,:2]);stat=np.full(n,np.nan);count=np.zeros(n,int)
                geom=(abs(targets[:,0])<180-r)&(abs(targets[:,1])<225-r)&(abs(targets[:,2]-offset)<depth/2)
                for i in np.flatnonzero(geom):
                    q=ix[tr.query_ball_point(targets[i,:2],r)];q=q[q!=ti[i]];count[i]=len(q)
                    if len(q)>=8:stat[i]=contrast2(coords[q,:2])
                key=f'a{angle}_o{offset}_d{depth}';arrays['shape_'+key]=stat;arrays['count_'+key]=count
                for name,cohort in cohorts.items():
                    valid=cohort&np.isfinite(stat);ll=valid&(stat>=cut)
                    views.append({'angle':angle,'offset':offset,'depth':depth,'cohort':name,'cohort_n':int(cohort.sum()),'inside_margin':int((cohort&geom).sum()),'valid':int(valid.sum()),'line_like':int(ll.sum()),'low_Pfil_line_like':int((ll&low).sum())})
        print('views complete',angle,flush=True)
    common_rows=[]
    for angle in [0,30,60,90]:
      for offset in [-5,0,5]:
       common=np.logical_and.reduce([np.isfinite(arrays[f'shape_a{angle}_o{offset}_d{dep}']) for dep in [5,10,20,40,80]])
       for name,cohort in cohorts.items():
        for dep in [5,10,20,40,80]:
         ll=cohort&common&(arrays[f'shape_a{angle}_o{offset}_d{dep}']>=cut)
         common_rows.append({'angle':angle,'offset':offset,'depth':dep,'cohort':name,'common_valid':int((cohort&common).sum()),'line_like':int(ll.sum())})
    args.out.mkdir(parents=True);np.savez_compressed(args.out/'geometry.npz',**arrays)
    result={'rows':rows,'views':views,'common_views':common_rows,'cohort_sizes':{k:int(v.sum()) for k,v in cohorts.items()},'job':os.environ['SLURM_JOB_ID'],
     'provenance':{'parent_sha256':digest(p),'audit_sha256':digest(args.audit),'script_sha256':digest(__file__),'geometry_sha256':digest(args.out/'geometry.npz')},
     'definitions':'L=(e1-e2)/e1, P=(e2-e3)/e1, S=e3/e1 for descending covariance eigenvalues; largest component defines descriptive line/sheet/isotropic dominance. Not tidal class probabilities. Self excluded by row index. Full/slab same 3D ball; slab abs(height)<10. Count-matched mean32 random subsets of full ball with slab count; random seed17017. Conditional 512 spatial-block bootstrap draws,30Mpc Cartesian blocks. Observed redshift-space tracers only; no correction for missing galaxies or peculiar velocities.'}
    (args.out/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('COMPLETE',args.out,flush=True)
if __name__=='__main__':main()
