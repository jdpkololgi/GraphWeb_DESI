"""Bounded old-sample neighbour check and training-field label parity."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import healpy as hp
from scipy.spatial import cKDTree
from astropy.cosmology import Planck18
from legacy_selection_tests_20261014 import OLD, OUT, ROOT, valid, sky
from legacy_label_scale_20261014 import member, CWEB
from desi_density_properties_20261014 import estimate, quartiles, dump
from desi_raw_environment_plots_20261014 import tied_quartiles

def main():
    torch.set_num_threads(1)
    archive=Path('/pscratch/sd/d/dkololgi/tng_illustris/cache/from_repo_root')
    field=member(CWEB,'cweb'); grid=np.memmap(CWEB,mode='r',dtype=field['dtype'],offset=field['offset'],shape=field['shape'])
    parity={}
    for name in ['processed_gcn_data_alpha_mc1e+09.pt','processed_gcn_data_alpha_scaler_saved_mc1e+09.pt']:
        data,features,targets=torch.load(archive/name,map_location='cpu',weights_only=False)
        pos=data.pos.numpy(); pixel=(pos/(205/.6774)*512).astype(int)
        labels=grid[pixel[:,0],pixel[:,1],pixel[:,2]]
        parity[name]=dict(rows=len(pos),field_label_mismatches=int(np.sum(labels!=data.y.numpy())),target_mismatches=int(np.sum(np.asarray(targets)!=data.y.numpy())),features=list(features.columns))
    dump(OUT/'training_label_parity.json',parity)
    old=dict(np.load(OLD/'legacy_cigale.npz')); good,mass,ssfr=valid(old)
    vac=pd.read_pickle(ROOT/'DESI_BGS_PRERELEASE_VAC.pkl');vac=vac[~vac.TARGETID.duplicated(False)]
    zg=np.linspace(0,.1,5001);dg=Planck18.comoving_distance(zg).value*Planck18.h
    def xyz(ra,dec,z):
        rr=np.interp(z,zg,dg); a=np.deg2rad(ra);d=np.deg2rad(dec)
        return rr[:,None]*np.stack([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)],1),rr
    xx,_=xyz(vac.RA.to_numpy(),vac.DEC.to_numpy(),vac.Z.to_numpy());tree=cKDTree(xx)
    q,rr=xyz(old['RA'],old['DEC'],old['z'])
    # Conservative rectangular-wedge and radial buffer; this is not a survey-mask correction.
    angles=np.deg2rad(np.stack([old['RA']-120,160-old['RA'],old['DEC']-14.5,30.6-old['DEC']],1))
    angles[:,:2]*=np.cos(np.deg2rad(old['DEC']))[:,None]
    boundary=np.minimum(np.min(rr[:,None]*np.sin(angles),axis=1),np.minimum(rr-np.interp(.01,zg,dg),np.interp(.06,zg,dg)-rr))
    ix=np.flatnonzero(good & sky(old['RA'],old['DEC']) & (boundary>7))
    counts=np.stack([tree.query_ball_point(q[ix],r,return_length=True,workers=1)-1 for r in [1,2,4,7]],1)
    # Independent brute force on 12 galaxies, all radii.
    for j in np.linspace(0,len(ix)-1,12).astype(int):
        dist=np.linalg.norm(xx-q[ix[j]],axis=1)
        assert np.array_equal(counts[j],[(dist<=r).sum()-1 for r in [1,2,4,7]])
    p=old['P_WEB'][ix];zb=np.floor(old['z'][ix]/.01).astype(int);reg=hp.ang2pix(8,old['RA'][ix],old['DEC'][ix],lonlat=True)
    weights={'soft_web':p}
    for j,r in enumerate([1,2,4,7]):weights['N'+str(r)],_=tied_quartiles(counts[:,j],zb)
    report=dict(training_parity=parity,query_rows=len(ix),tracers=len(vac),median_counts=np.median(counts,axis=0).tolist(),zero_fraction={str(r):float((counts[:,j]==0).mean()) for j,r in enumerate([1,2,4,7])},methods={},limits='Same old galaxies across estimators; deduplicated saved mass-cut catalogue as tracers. Rectangular wedge and radial 7 Mpc/h buffer, no survey angular-completeness correction. Density quartiles use tie-equivalent fractional memberships, matching the current P12-A count analysis. They are different groups from web classes. No causal or model-accuracy claim.')
    for width in [.1,.2]:
        cells=np.floor((mass[ix]-9)/width).astype(int)*100+zb
        try:
            result,_=estimate((ssfr[ix]<-11).astype(float),weights,cells,reg,np.zeros(len(ix)),minimum=20,nboot=64)
            report['methods'][str(width)]=result
        except (ValueError,AssertionError) as e:report['methods'][str(width)]={'insufficient_support':str(e)}
    dump(OUT/'old_density_scale.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
