"""Read-only full-table selection ablations; no posterior fitting or new phase access."""
from pathlib import Path
import sys,json
import numpy as np
import fitsio,healpy as hp
import loa_trial as t

ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
MOCK='/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl6/kibo-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.fits'

def scan(path,kind,common):
    counts={};sums={};total=0;zwarn={};photo={};maghist=np.zeros((8,100),dtype='i8');edges=np.linspace(15,21,101)
    def add(name,mask,binid):
        counts.setdefault(name,np.zeros(8,dtype='i8'));counts[name]+=np.bincount(binid[mask],minlength=8)[:8]
    with fitsio.FITS(path) as f:
      for start in range(0,f[1].get_nrows(),250000):
        d=f[1][start:min(start+250000,f[1].get_nrows())];total+=len(d)
        if kind=='mock':
            tz=d['TRUEZ'];ts=np.searchsorted([.15,.25,.35,.45,.55],tz,side='right')-1;tv=np.isfinite(tz)&(ts>=0)&(ts<4)
            tc=t.galactic_cap(d['RA'],d['DEC']);tb=np.clip(tc*4+ts,0,7)
            add('TRUEZ_all_targets',tv,tb);add('TRUEZ_zwarn0',tv&(d['ZWARN']==0),tb);add('TRUEZ_assigned',tv&d['LOCATION_ASSIGNED'].astype(bool),tb)
        z=d['Z'] if kind=='observed_mock' else d['Z_not4clus'];s=np.searchsorted([.15,.25,.35,.45,.55],z,side='right')-1;valid=np.isfinite(z)&(s>=0)&(s<4);d=d[valid];z=z[valid];s=s[valid]
        cap=t.galactic_cap(d['RA'],d['DEC']);bins=cap*4+s;pix=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True);zgood=d['ZWARN']==0
        masks={'redshift_range':np.ones(len(d),bool),'zwarn0':zgood}
        if kind=='loa':
            galaxy=np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY';base=zgood&galaxy&(d['DELTACHI2']>=25)
            masks.update(base=base,zwarn0_galaxy=zgood&galaxy,base_dchi40=base&(d['DELTACHI2']>40),lss40=zgood&(d['DELTACHI2']>40),base_dchi100=base&(d['DELTACHI2']>100),base_dchi1000=base&(d['DELTACHI2']>1000))
            flux=d['FLUX_R']/d['MW_TRANSMISSION_R'];r=np.full(len(d),np.nan);ok=flux>0;r[ok]=22.5-2.5*np.log10(flux[ok])
            masks.update(base_r19p5=base&(r<19.5),base_r19p54=base&(r<19.54),base_fiberstatus0=base&(d['COADD_FIBERSTATUS']==0))
            for field in ['WEIGHT_ZFAIL','mod_success_rate','PROB_OBS']:
                v=d[field];ok=base&np.isfinite(v)&(v>0);sums.setdefault(field,dict(sum=np.zeros(8),rows=np.zeros(8,dtype='i8')));sums[field]['sum']+=np.bincount(bins[ok],weights=v[ok],minlength=8);sums[field]['rows']+=np.bincount(bins[ok],minlength=8)
        else:
            base=zgood;masks['base']=base
            if kind=='observed_mock':
                r=d['R_MAG_APP'];masks.update(base_r19p5=base&(r<19.5),base_r19p54=base&(r<19.54),base_valid_box=base&(d['BOX_INDEX']>=0))
            else:r=np.full(len(d),np.nan)
        for name in ['GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED']:
            masks['base_'+name]=base&d[name].astype(bool)
        masks.update(base_bright_bit=base&((d['BGS_TARGET']&2)!=0),base_common_pixels=base&common[pix],base_maskbits0=base&(d['MASKBITS']==0),base_nobs=base&(d['NOBS_G']>0)&(d['NOBS_R']>0)&(d['NOBS_Z']>0))
        if kind=='loa':
            masks['common_r19p5_dchi40']=base&common[pix]&(r<19.5)&(d['DELTACHI2']>40)
            masks['base_N_rge19p5']=base&(np.char.strip(d['PHOTSYS'].astype('U'))=='N')&(r>=19.5)
            masks['base_S_rge19p5']=base&(np.char.strip(d['PHOTSYS'].astype('U'))=='S')&(r>=19.5)
        for n in [1,2,3]:masks[f'base_ntile_ge{n}']=base&(d['NTILE']>=n)
        for label,mask in masks.items():add(label,mask,bins)
        for c in [b'N',b'S','N','S']:
            sel=base&(d['PHOTSYS']==c)
            if sel.any():add('base_photsys_'+str(c),sel,bins)
        for b in range(8):maghist[b]+=np.histogram(r[base&(bins==b)],edges)[0]
        vals,num=np.unique(d['ZWARN'],return_counts=True)
        for v,n in zip(vals,num):zwarn[str(v)]=zwarn.get(str(v),0)+int(n)
        if kind=='mock':
            # TRUEZ counts include unobserved targets where the supplied column is valid.
            tz=d['TRUEZ'];truthvalid=np.isfinite(tz)&(tz>=.15)&(tz<.55)
            add('base_TRUEZ_matches_Z',base&truthvalid&(np.abs(tz-z)<.02),bins)
    return dict(path=path,total_rows=total,count_order='SGC shells0-3 then NGC shells0-3',counts={k:v.tolist() for k,v in counts.items()},weight_summary={k:{kk:vv.tolist() for kk,vv in v.items()} for k,v in sums.items()},zwarn_in_science_z=zwarn,rmag_edges=edges.tolist(),rmag_hist=maghist.tolist())

def main():
    loa=t.read(ROOT/'INPUTS_READY.json')['source']['path'];gal=np.load(ROOT/'galaxy_angular_support.npy');mock=np.load(t.B/'ph006/p3_fields/angular_support_nside256.npz')['support'];common=gal&mock
    paths=[('loa',loa),('mock',MOCK),('observed_mock',str(t.B/'ph006/catalogues/observed/ph006_bgs_bright_full_observed_with_tweb.fits'))];out={}
    for kind,path in paths:
        out[kind]=scan(path,kind,common);print(kind,json.dumps(out[kind]['counts']),flush=True)
    t.save(ROOT/'SELECTION_AUDIT_V2.json',dict(catalogues=out,common_pixels=int(common.sum()),note='Mock has no DELTACHI2/SPECTYPE. Absent/noisy-redshift models cannot be made equivalent by applying unavailable columns. Individual ablations are diagnostic, not newly chosen science cuts.'))
if __name__=='__main__':main()
