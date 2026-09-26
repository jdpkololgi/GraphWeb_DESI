"""Bounded ph000 producer schema and colour-mapping probe; no count inference."""
import hashlib
import json
from pathlib import Path
import fitsio
import h5py
import numpy as np
from photometry_recipe_fingerprint import k, stats, OUT

ROOT = Path('/global/cfs/cdirs/desi/users/jpiat/abacus_mocks/cutsky/z0.200/AbacusSummit_base_c000/ph000')
DEST = Path(__file__).resolve().parents[2] / 'docs/evidence/p12a_alignment_execution_20260926/PRODUCER_PROBE.json'

def main():
    result = {'scope': 'Only exposed ph000; 2048 evenly spaced rows per file. Not a population census.', 'files': [], 'holi': []}
    for p in sorted(ROOT.glob('*.fits')) + sorted((ROOT/'forFA').glob('*.fits')):
        with fitsio.FITS(p) as f:
            t = f[1]; n = t.get_nrows(); names = t.get_colnames()
            rows = np.linspace(0, n-1, 2048, dtype=np.int64)
            cols = ['RA','DEC','R_MAG_APP','G_R_OBS','G_R_REST']
            zn = 'Z' if 'Z' in names else 'RSDZ'
            cols += [zn] + [s for s in ['BGS_TARGET','TARGETID','ZWARN','PRIORITY'] if s in names]
            d = t.read(rows=rows, columns=cols)
            z=d[zn]; c=d['G_R_REST']; keep=(z>=.15)&(z<.55)&(d['R_MAG_APP']<19.5)
            rec = dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,rows=n,columns=names,
                       primary_header=str(f[0].read_header()),table_header=str(t.read_header()),
                       sample_sha256=hashlib.sha256(d.tobytes()).hexdigest(),selected_sample_rows=int(keep.sum()),
                       sample_ra_range=[float(d['RA'].min()),float(d['RA'].max())],
                       sample_dec_range=[float(d['DEC'].min()),float(d['DEC'].max())],
                       old_colour_mapping_residual=stats((d['G_R_OBS']-c-k(z,c,'g')+k(z,c,'r'))[keep]))
            for field in ['BGS_TARGET','ZWARN','PRIORITY']:
                if field in names:
                    u,counts=np.unique(d[field],return_counts=True)
                    rec[field+'_sample_counts']={str(int(a)):int(b) for a,b in zip(u,counts)}
            result['files'].append(rec)
    base=Path('/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.82/seed0000')
    for name in ['holi_BGS_v4.82_GCcomb_clustering.dat.h5','holi_BGS-NONKP_v4.82_GCcomb_clustering.dat.h5']:
        p=base/name
        with h5py.File(p) as h:
            result['holi'].append(dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,
                fields={key:dict(shape=list(h[key].shape),dtype=str(h[key].dtype)) for key in h},
                root_attributes={key:str(val) for key,val in h.attrs.items()}))
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['table_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (OUT/'source').glob('k_corr*')}
    DEST.write_text(json.dumps(result,indent=2)+'\n')
    for r in result['files']:
        print(Path(r['path']).name,r['rows'],r['selected_sample_rows'],r['old_colour_mapping_residual'],r.get('BGS_TARGET_sample_counts'))

if __name__ == '__main__':
    main()
