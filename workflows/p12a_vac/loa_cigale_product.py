"""Carry the validated CIGALE join onto every original posterior-VAC column."""
import argparse,json
from pathlib import Path
import numpy as np
import fitsio
from property_environment import VAC,digest

def main(root):
    receipt=json.loads((root/'JOIN.json').read_text());p=Path(receipt['joined_path'])
    assert digest(p)==receipt['joined_sha256'] and digest(VAC)==receipt['vac_sha256']
    c=fitsio.read(p);d=fitsio.read(VAC);assert np.array_equal(c['TARGETID'],d['TARGETID'])
    extra=[n for n in c.dtype.names if n not in d.dtype.names]
    output=np.empty(len(d),dtype=list(d.dtype.descr)+[(n,c.dtype[n]) for n in extra])
    for n in d.dtype.names:output[n]=d[n]
    for n in extra:output[n]=c[n]
    with fitsio.FITS(VAC) as f:
        header=f[1].read_header();original_units={n:header.get(f'TUNIT{i+1}','') for i,n in enumerate(d.dtype.names)}
        (root/'ORIGINAL_VAC_HEADER.txt').write_text(str(header))
    units=[]
    for n in output.dtype.names:
        unit=original_units.get(n,'')
        if '_CG_' in n:
            if n.startswith(('MASS','MASSERR')):unit='solMass'
            elif n.startswith(('SFR','SFRERR')):unit='solMass yr-1'
            elif n.startswith(('AGE','AGEERR')):unit='Myr'
            elif n.startswith(('AV','AVERR')):unit='mag'
        if n=='CIGALE_SEP_ARCSEC':unit='arcsec'
        units.append(unit)
    dest=root/'DESI_LOA_P12A_CIGALE_FULL_VAC.fits'
    fitsio.write(dest,output,clobber=True,units=units,header={'PROVIS':True,'VACSHA':receipt['vac_sha256'],'CGSHA':receipt['source_sha256'],'CGJOIN':'TARGETID, sep<1 arcsec, abs(dz)/(1+z)<0.001; duplicates excluded'})
    check=fitsio.read(dest,columns=list(d.dtype.names))
    for n in d.dtype.names:
        assert np.array_equal(d[n],check[n],equal_nan=d[n].dtype.kind in 'fc'),n
    r={'path':str(dest),'sha256':digest(dest),'rows':len(d),'original_columns':list(d.dtype.names),'cigale_columns':extra,'all_original_columns_exactly_preserved':True,'join_receipt_sha256':digest(root/'JOIN.json'),'script_sha256':digest(__file__)}
    (root/'PRODUCT.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
