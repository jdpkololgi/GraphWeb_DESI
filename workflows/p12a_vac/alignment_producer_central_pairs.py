"""Bounded shared-central-location comparison, not proof of galaxy identity."""
from pathlib import Path
import hashlib,json
import fitsio
import numpy as np
from alignment_producer_probe import ROOT,DEST
from photometry_recipe_fingerprint import stats

cols=['HALO_ID','CEN','RA','DEC','Z','R_MAG_APP','R_MAG_ABS','G_R_OBS','G_R_REST']
a={}
for branch in ['N','S']:
    with fitsio.FITS(ROOT/f'galaxy_cut_sky_{branch}.fits') as f:
        d=f[1].read(rows=np.arange(10000),columns=cols)
    d=d[d['CEN']==1]
    keys=[(int(t['HALO_ID']),float(t['RA']),float(t['DEC']),float(t['Z'])) for t in d]
    assert len(keys)==len(set(keys)), 'Ambiguous central-location keys'
    a[branch]=dict(zip(keys,d))
keys=sorted(set(a['N'])&set(a['S']))
res={'contract':'First 10000 raw rows per producer branch, CEN==1, exact HALO_ID/RA/DEC/Z key; shared central locations, not verified galaxy IDs. Not a representative population sample.',
     'centrals':{b:len(a[b]) for b in a},'matched':len(keys),'deltas_S_minus_N':{},
     'source_stat':{b:[(ROOT/f'galaxy_cut_sky_{b}.fits').stat().st_size,(ROOT/f'galaxy_cut_sky_{b}.fits').stat().st_mtime_ns] for b in a},
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
for c in ['R_MAG_APP','R_MAG_ABS','G_R_OBS','G_R_REST']:
    res['deltas_S_minus_N'][c]=stats(np.array([float(a['S'][key][c])-float(a['N'][key][c]) for key in keys]))
(DEST.parent/'PRODUCER_CENTRAL_PAIRS.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
