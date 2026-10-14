#!/usr/bin/env python3
"""New non-neural Loa wedge environments; single-thread, bounded neighbour queries."""
from pathlib import Path
import argparse,json,hashlib,time,os
import numpy as np
import fitsio,healpy as hp
from astropy.cosmology import Planck18
from scipy.spatial import cKDTree
from numba import njit
ROOT=Path('/global/homes/d/dkololgi/GraphWeb_DESI')
VAC=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1/DESI_LOA_P12A_HALO48_FULL_SURVEY_VAC.fits')
INPUT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
PRODUCT=ROOT/'docs/figures/p12a_loa_cigale_20260929/PRODUCT.json'
H=.6766;LIGHT=299792.458;VMAX=1000.;RPMAX=10/H;RADII=np.array([1.,2.,4.,7.])/H

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def unit(ra,dec):
 a=np.deg2rad(ra);b=np.deg2rad(dec);return np.column_stack((np.cos(a)*np.cos(b),np.sin(a)*np.cos(b),np.sin(b)))

@njit(cache=True)
def reduce_neighbors(qidx,ptr,nei,xyz,z,ids,radius,u):
 n=len(qidx);counts=np.zeros((n,4),np.int64);out=np.full((n,7),np.nan)
 # output d5, projected_count, Q3_unbiased, Qsky_unbiased, LOS_P2, Ndir3, Ndir2
 for a in range(n):
  i=qidx[a];near5=np.full(5,np.inf);npv=0;nn=0;nt=0;m=np.zeros((3,3));c2=0.;s2=0.;los=0.
  east=np.array([-u[i,1],u[i,0],0.]);east/=np.sqrt((east*east).sum());north=np.cross(u[i],east)
  for pp in range(ptr[a],ptr[a+1]):
   j=nei[pp]
   if ids[j]==ids[i]:continue
   v=xyz[j]-xyz[i];r2=(v*v).sum()
   for k in range(4):
    if r2<=RADII[k]**2:counts[a,k]+=1
   if r2<=RADII[3]**2 and r2>1e-20:
    vv=v/np.sqrt(r2);nn+=1
    for b in range(3):
     for c in range(3):m[b,c]+=vv[b]*vv[c]
    mu=(vv*u[i]).sum();los+=(3*mu*mu-1)/2
    vx=(v*east).sum();vy=(v*north).sum();vt=vx*vx+vy*vy
    if vt>1e-20:nt+=1;c2+=(vx*vx-vy*vy)/vt;s2+=2*vx*vy/vt
   if LIGHT*abs(z[j]-z[i])/(1+z[i])<=VMAX:
    chord=np.sqrt(((u[j]-u[i])**2).sum());theta=2*np.arcsin(min(1.,chord/2));rp=radius[i]*theta
    if rp<=RPMAX:
     npv+=1
     if rp<near5[4]:
      near5[4]=rp
      for k in range(4,0,-1):
       if near5[k]<near5[k-1]:tmp=near5[k];near5[k]=near5[k-1];near5[k-1]=tmp
  if npv>=5:out[a,0]=near5[4]
  out[a,1]=npv;out[a,5]=nn;out[a,6]=nt
  if nn>1:
   tr=(m*m).sum()/nn**2;raw=(3*tr-1)/2;out[a,2]=(nn*raw-1)/(nn-1);out[a,4]=los/nn
  if nt>1:
   raw=(c2*c2+s2*s2)/nt**2;out[a,3]=(nt*raw-1)/(nt-1)
 return counts,out

def build(out):
 t=time.monotonic();receipt=json.loads(PRODUCT.read_text());ready=json.loads((INPUT/'INPUTS_READY.json').read_text())
 for name in ['random_angular.npz','boundary_angles.npz']:assert sha(INPUT/name)==ready['arrays_sha256'][name]
 cols=['TARGETID','RA','DEC','Z','CAP','SUPPORTED','BOUNDARY_MPC','NTILDE_MPC3']
 d=fitsio.read(VAC,columns=cols);assert len(d)==receipt['rows'];print('Read',len(d),'Loa tracers',flush=True)
 q=(d['RA']>=120)&(d['RA']<160)&(d['DEC']>=14.5)&(d['DEC']<30.6)&(d['Z']>=.2)&(d['Z']<.3)&d['SUPPORTED'];rows=np.flatnonzero(q);assert len(rows)>100000
 zg=np.linspace(.1,.6,10001);rg=Planck18.comoving_distance(zg).value
 r=np.interp(d['Z'],zg,rg);uu=unit(d['RA'],d['DEC']);xyz=uu*r[:,None]
 zq=d['Z'][rows];rq=r[rows];dz=VMAX/LIGHT*(1+zq);rhi=np.interp(zq+dz,zg,rg);rlo=np.interp(zq-dz,zg,rg)
 # Euclidean ball guaranteed to enclose the projected aperture and velocity window.
 search=np.sqrt(np.maximum(rhi-rq,rq-rlo)**2+(rhi/rq)*RPMAX**2)+1e-5
 pad=max(search.max(),RADII[-1])+1e-3;lo=xyz[rows].min(0)-pad;hi=xyz[rows].max(0)+pad
 take=np.all((xyz>=lo)&(xyz<=hi),axis=1);tracer_rows=np.flatnonzero(take);inv=np.full(len(d),-1,dtype=np.int64);inv[tracer_rows]=np.arange(len(tracer_rows));qidx=inv[rows];assert (qidx>=0).all()
 dd=d[rows];tr=d[take];xx=xyz[take];ut=uu[take];rt=r[take];del d,xyz,uu,r,inv
 tz=tr['Z'].astype(np.float64);tid=tr['TARGETID'].astype(np.int64)
 tree=cKDTree(xx);counts=np.zeros((len(rows),4),dtype=np.int64);metrics=np.zeros((len(rows),7));candidates=0
 for start in range(0,len(rows),256):
  stop=min(start+256,len(rows));ii=qidx[start:stop]
  lists=tree.query_ball_point(xx[ii],search[start:stop],workers=1,return_sorted=False)
  lens=np.array([len(x) for x in lists]);ptr=np.r_[0,np.cumsum(lens)];nei=np.concatenate(lists).astype(np.int64);candidates+=len(nei)
  counts[start:stop],metrics[start:stop]=reduce_neighbors(ii,ptr,nei,xx,tz,tid,rt,ut)
  if start%10240==0:print('Measured',stop,'/',len(rows),'elapsed',round(time.monotonic()-t,1),flush=True)
 # Independent brute-force checks on 12 deterministic galaxies, including projected search completeness.
 rng=np.random.default_rng(20261014);tests=np.sort(rng.choice(len(rows),12,replace=False));errors=[]
 for a in tests:
  i=qidx[a];v=xx-xx[i];r2=(v*v).sum(1);other=tr['TARGETID']!=tr['TARGETID'][i]
  expected=np.array([np.sum(other&(r2<=rr**2)) for rr in RADII]);assert np.array_equal(expected,counts[a])
  theta=2*np.arcsin(np.clip(np.linalg.norm(ut-ut[i],axis=1)/2,0,1));rp=rt[i]*theta;k=other&(LIGHT*np.abs(tr['Z']-tr['Z'][i])/(1+tr['Z'][i])<=VMAX)&(rp<=RPMAX)
  assert k.sum()==metrics[a,1]
  if k.sum()>=5:assert np.isclose(np.partition(rp[k],4)[4],metrics[a,0],rtol=1e-10,atol=1e-10)
  errors.append(int(a))
 # Exact ID parity to the augmented current catalogue; no historical property-cache selection.
 ccols=['TARGETID','Z','P_WEB','EIGENVALUE_MEAN','CIGALE_MATCHED','CIGALE_SPECTYPE','MASS_CG_15','SFR_CG_15','AGE_CG_15','AV_CG_15','FAGN_CG_15']
 cg=fitsio.read(receipt['path'],rows=rows,columns=ccols);assert np.array_equal(cg['TARGETID'],dd['TARGETID']) and np.array_equal(cg['Z'],dd['Z'])
 angular=np.load(INPUT/'random_angular.npz');ang=np.load(INPUT/'boundary_angles.npz');pix=hp.ang2pix(256,dd['RA'],dd['DEC'],lonlat=True)
 angle=np.where(dd['CAP']==1,ang['NGC'][pix],ang['SGC'][pix]);lower=np.maximum(angle-2*hp.max_pixrad(256),0);boundary=rq*np.sin(np.minimum(lower,np.pi/2))
 supported=angular['support'][pix]&((angular['domain'][pix]//2)==dd['CAP'])
 # The VAC tracer radial limits .15-.55 comfortably enclose all these apertures.
 radedge=np.minimum(rq-np.interp(.15,zg,rg),np.interp(.55,zg,rg)-rq)
 valid_sphere=supported&(boundary>RADII[-1])&(radedge>RADII[-1]);valid_proj=supported&(lower>metrics[:,0]/rq)&(radedge>np.maximum(rhi-rq,rq-rlo))&np.isfinite(metrics[:,0])&(metrics[:,0]>0)
 ratio=np.divide(counts[:,1]*(RADII[3]**3-RADII[1]**3), (counts[:,3]-counts[:,1])*RADII[1]**3,out=np.full(len(rows),np.nan),where=counts[:,3]>counts[:,1])
 sigma=5/(np.pi*(metrics[:,0]*H)**2)
 np.savez_compressed(out/'raw_environment_catalogue.npz',catalogue_rows=rows,TARGETID=dd['TARGETID'],RA=dd['RA'],DEC=dd['DEC'],z=dd['Z'],cap=dd['CAP'],NTILDE_MPC3=dd['NTILDE_MPC3'],counts=counts,radii_hMpc=RADII*H,d5_hMpc=metrics[:,0]*H,sigma5=sigma,projected_candidates=metrics[:,1].astype(int),Q3=metrics[:,2],Qsky=metrics[:,3],LOS_P2=metrics[:,4],Ndir3=metrics[:,5].astype(int),Ndir2=metrics[:,6].astype(int),inner_outer_ratio=ratio,boundary_lower_Mpc=boundary,valid_sphere=valid_sphere,valid_sigma5=valid_proj,region=hp.ang2pix(8,dd['RA'],dd['DEC'],lonlat=True),**{n:cg[n] for n in cg.dtype.names if n not in ['TARGETID','Z']})
 dump(out/'measurement_provenance.json',{'source':str(VAC),'source_receipt':receipt,'full_fits_hash_recomputed':False,'source_binding':'Frozen receipts, row count, exact TARGETID/redshift crosswalk; fresh extracted-output hash.','support_hashes':{n:sha(INPUT/n) for n in ['INPUTS_READY.json','random_angular.npz','boundary_angles.npz']},'script_sha256':sha(__file__),'catalogue_sha256':sha(out/'raw_environment_catalogue.npz'),'queries':len(rows),'padded_tracers':len(tr),'candidate_pairs_processed':candidates,'elapsed_seconds':time.monotonic()-t,'selection':'SUPPORTED Loa RA[120,160),DEC[14.5,30.6),z[.20,.30); all selected rows measured before CIGALE cuts. Full Loa VAC neighbours within padded Cartesian bounding box, including outside wedge.','coordinates':'Planck18 comoving Mpc using 10001-point interpolation z=.1-.6; h=.6766.','sigma5':'5/(pi*d5^2), d5 in comoving h^-1 Mpc at target distance, angular great-circle separation; |c dz/(1+z_target)| <=1000 km/s; self excluded by TARGETID; search limited to10 h^-1Mpc, missing fifth neighbour flagged NaN.','counts':'Self-excluded redshift-space spheres R=1,2,4,7 h^-1Mpc. Inner/outer density ratio compares R<2 with2<R<7; no pseudocount.','Q3':'M=mean(unit-neighbour-vector outer product) within R7; A=(3tr(M^2)-1)/2; Q3=(N*A-1)/(N-1). Negative debiased values retained. A tests quadrupole anisotropy only.','Qsky':'Same R7 neighbours projected onto target sky plane; A2=|mean(exp(2i phi))|^2; Qsky=(N*A2-1)/(N-1). LOS_P2=mean(P2(unit_neighbour dot line_of_sight)).','shape_validity':'Statistics recorded for Ndir>1; property plots require Ndir>=20. Finite-N subtraction assumes independent isotropic directions; clustered correlations and redshift distortions remain.','boundaries':'Use stored nside256 angular distance at actual galaxy HEALPix, subtract2*max_pixrad for conservative finite-pixel buffer; exact target radius, mask/cap support. For spheres require lower distance>R7; projected require lower angle>d5/chi and full velocity window inside tracer radial limits. This addresses geometry, not variable angular completeness.','qa':{'counts_and_sigma5_bruteforce_queries':errors,'ids_unique':bool(len(np.unique(dd['TARGETID']))==len(dd)),'count_monotonic':bool(np.all(np.diff(counts,axis=1)>=0)),'valid_sphere_rows':int(valid_sphere.sum()),'valid_sigma5_rows':int(valid_proj.sum()),'fewer_than_five_projected':int(np.sum(~np.isfinite(metrics[:,0])))},'execution':'Single-thread login execution requested by user; bounded wedge only, no Slurm changes or ML inference.'})
 print('COMPLETE',len(rows),'queries;',len(tr),'tracers;',round(time.monotonic()-t,1),'seconds',flush=True)

def test():
 # Exact self exclusion, known fifth projected radius, finite-N formula, rotation invariance.
 xyz=np.array([[100.,0,0],[100,1,0],[100,-1,0],[100,0,1],[100,0,-1],[101,0,0],[99,0,0]],float);r=np.linalg.norm(xyz,axis=1);u=xyz/r[:,None];z=np.full(7,.025);ids=np.arange(7,dtype=np.int64);q=np.array([0]);ptr=np.array([0,7]);nei=np.arange(7)
 c,m=reduce_neighbors(q,ptr,nei,xyz,z,ids,r,u);assert np.all(c==6);assert np.isclose(m[0,2],-.2)
 assert np.isclose(m[0,0],100*np.arctan(.01))
 rot=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
 c2,m2=reduce_neighbors(q,ptr,nei,xyz@rot,z,ids,r,u@rot);assert np.array_equal(c,c2) and np.allclose(m,m2,equal_nan=True)
 # Collinear directions yield Q3=1; isotropic iid ensemble has zero expected debiased Q3.
 rng=np.random.default_rng(42);v=rng.normal(size=(2000,20,3));v/=np.linalg.norm(v,axis=2,keepdims=True);M=np.einsum('bni,bnj->bij',v,v)/20;A=(3*(M*M).sum((1,2))-1)/2;assert abs(np.mean((20*A-1)/19))<.003
 print('Self exclusion, known projected d5, finite-N bias, and rotation tests passed')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path);p.add_argument('--test',action='store_true');a=p.parse_args()
 if a.test:test()
 else:a.out.mkdir(parents=True,exist_ok=True);build(a.out)
