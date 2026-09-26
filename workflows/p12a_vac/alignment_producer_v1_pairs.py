"""Bounded canonical-v1 central-location crosswalk to producer N/S."""
import hashlib,json
from pathlib import Path
import fitsio
import numpy as np
from alignment_producer_probe import ROOT,DEST
from alignment_parent_screen import BASE
from photometry_recipe_fingerprint import stats

paths={b:ROOT/f'galaxy_cut_sky_{b}.fits' for b in ['N','S']}
paths['v1']=BASE/'v1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
fields=['R_MAG_APP','R_MAG_ABS','G_R_OBS','G_R_REST']
a={}
for b,p in paths.items():
    with fitsio.FITS(p) as f:d=f[1].read(rows=np.arange(10000),columns=['HALO_ID','CEN','RA','DEC','Z','Z_COSMO','vx','vy','vz']+fields)
    d=d[d['CEN']==1]
    keys=[(int(t['HALO_ID']),float(t['RA']),float(t['DEC']),float(t['Z'])) for t in d]
    assert len(keys)==len(set(keys))
    a[b]=dict(zip(keys,d))
r={'contract':'First 10000 raw rows, central-location exact HALO_ID/RA/DEC/Z key; prefix-biased sample, not full identity or generation provenance.', 'pairs':{},
'source_stat':{b:[p.stat().st_size,p.stat().st_mtime_ns] for b,p in paths.items()},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
for b in ['N','S']:
    keys=sorted(set(a[b])&set(a['v1']))
    r['pairs'][b]={'matched':len(keys),'deltas_producer_minus_v1':{},'all_four_photometry_fields_exact':sum(all(a[b][key][f]==a['v1'][key][f] for f in fields) for key in keys)}
    for f in fields:r['pairs'][b]['deltas_producer_minus_v1'][f]=stats(np.array([float(a[b][key][f])-float(a['v1'][key][f]) for key in keys]))
r['angular_only_pairs']={}
bb={b:{key[:-1]:value for key,value in rows.items()} for b,rows in a.items()}
for b in a:assert len(bb[b])==len(a[b])
for b in ['N','S']:
    keys=sorted(set(bb[b])&set(bb['v1']))
    r['angular_only_pairs'][b]={'matched':len(keys),'deltas_producer_minus_v1':{}}
    for f in fields+['Z','Z_COSMO','vx','vy','vz']:
        r['angular_only_pairs'][b]['deltas_producer_minus_v1'][f]=stats(np.array([float(bb[b][key][f])-float(bb['v1'][key][f]) for key in keys]))
(DEST.parent/'PRODUCER_V1_PAIRS.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
