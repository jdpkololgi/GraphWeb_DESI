"""Crosswalk paired Uchuu Loa clustering selection to our actual full sample."""
import argparse
import os
from pathlib import Path
import fitsio
import healpy as hp
import numpy as np
from alignment_uchuu_dr2_screen import PATHS
from flag_joint_magnitude import REAL
from alignment_parent_screen import ZE
from mock_selection_test import cap,digest,save

def run(root,parent):
    if (root/'COMPLETE.json').exists():raise FileExistsError(root)
    mask=np.load(parent/'common.npy')
    d=fitsio.read(PATHS['paired_loa'],ext=1,columns=['TARGETID','Z','RA','DEC'])
    order=np.argsort(d['TARGETID']);d=d[order];assert np.all(np.diff(d['TARGETID'])>0)
    ids=d['TARGETID'];pcap=cap(d['RA'],d['DEC']);pkeep=(d['Z']>=.15)&(d['Z']<.55)&mask[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
    seen=np.zeros(len(d),bool)
    counts={k:np.zeros((2,40),dtype='i8') for k in ['current_loa','paired_fail_zwarn','paired_fail_dchi25','paired_not_galaxy','paired_current_z_outside','paired_any_current_science_fail']}
    stats=dict(paired_matched_ids=0)
    before=[REAL.stat().st_size,REAL.stat().st_mtime_ns]
    with fitsio.FITS(REAL) as f:
        n=f[1].get_nrows()
        for start in range(0,n,250000):
            a=f[1].read(rows=np.arange(start,min(start+250000,n)),columns=['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE'])
            z=a['Z_not4clus'];c=cap(a['RA'],a['DEC']);sky=mask[hp.ang2pix(256,a['RA'],a['DEC'],lonlat=True)]
            zw=a['ZWARN']==0;dc=a['DELTACHI2']>=25;gal=np.char.strip(a['SPECTYPE'].astype('U'))=='GALAXY';zr=(z>=.15)&(z<.55);q=zw&dc&gal&zr
            for i in (0,1):counts['current_loa'][i]+=np.histogram(z[q&sky&(c==i)],bins=ZE)[0]
            ix=np.searchsorted(ids,a['TARGETID']);matched=ix<len(ids);ix=np.minimum(ix,len(ids)-1);matched &=ids[ix]==a['TARGETID']
            seen[ix[matched]]=True;stats['paired_matched_ids']+=int(matched.sum())
            for key,fail in [('paired_fail_zwarn',~zw),('paired_fail_dchi25',~dc),('paired_not_galaxy',~gal),('paired_current_z_outside',~zr),('paired_any_current_science_fail',~q)]:
                use=matched&pkeep[ix]&fail
                for i in (0,1):
                    kk=use&(pcap[ix]==i);counts[key][i]+=np.histogram(d['Z'][ix[kk]],bins=ZE)[0]
            if start%2000000==0:print('current Loa parity',start,flush=True)
    assert before==[REAL.stat().st_size,REAL.stat().st_mtime_ns]
    counts['paired_missing_current']=np.array([np.histogram(d['Z'][pkeep&~seen&(pcap==i)],bins=ZE)[0] for i in (0,1)])
    mock=np.load(parent/'HISTOGRAMS.npz')['uchuu_altmtl'];broad={k:v.reshape(2,4,10).sum(-1) for k,v in counts.items()}
    np.savez_compressed(root/'HISTOGRAMS.npz',**counts,z_edges=ZE,mock=mock)
    save(root/'COMPLETE.json',dict(counts={k:v.tolist() for k,v in broad.items()},stats=stats,
        mock_to_current_loa=(mock.reshape(2,4,10).sum(-1)/broad['current_loa']).tolist(),
        sources={str(REAL):before,str(PATHS['paired_loa']):[PATHS['paired_loa'].stat().st_size,PATHS['paired_loa'].stat().st_mtime_ns]},
        script_sha256=digest(__file__),histogram_sha256=digest(root/'HISTOGRAMS.npz'),job=os.environ.get('SLURM_JOB_ID'),
        notes=['Paired failures are nonexclusive; bins use paired Z to account for that sample. Current baseline counts use current Z.',
               'Current full data quality is ZWARN0, DELTACHI2>=25, SPECTYPE=GALAXY, .15<=Z_not4clus<.55.',
               'Mock processed clustering selection not rederived from parent; this improves the data-side comparison without proving full observation parity.']))
    print('COMPLETE parity',mock.reshape(2,4,10).sum(-1)/broad['current_loa'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--parent',type=Path,required=True)
    a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True);run(a.root,a.parent)
