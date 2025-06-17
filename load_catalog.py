import astropy.constants as c
import astropy.units as u
from astropy.table import Table, join, vstack, unique, hstack, unique
from astropy.coordinates import SkyCoord
from astropy.io import fits
import fitsio
import os
import numpy as np

## Set up the path to the fastscpec fit catalogs, and create the array with the list of all the various fastspec-loa catalogs 
## The list below is correct as of 28.3.25 for v1 of the Loa processing 
## This includes the patch to fix the emission line flux uncertaintites 

# fastspec_path = '/pscratch/sd/i/ioannis/fastspecfit/data/loa-bug/catalogs/'
fastspec_path = '/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/'

fastspec_catalogs = [
    'fastspec-loa-cmx-other.fits',
    'fastspec-loa-main-backup.fits',
    'fastspec-loa-main-bright-nside1-hp00.fits',
    'fastspec-loa-main-bright-nside1-hp01.fits',
    'fastspec-loa-main-bright-nside1-hp02.fits',
    'fastspec-loa-main-bright-nside1-hp03.fits',
    'fastspec-loa-main-bright-nside1-hp04.fits',
    'fastspec-loa-main-bright-nside1-hp05.fits',
    'fastspec-loa-main-bright-nside1-hp06.fits',
    'fastspec-loa-main-bright-nside1-hp07.fits',
    'fastspec-loa-main-bright-nside1-hp08.fits',
    'fastspec-loa-main-bright-nside1-hp09.fits',
    'fastspec-loa-main-bright-nside1-hp10.fits',
    'fastspec-loa-main-bright-nside1-hp11.fits',
    'fastspec-loa-main-dark-nside1-hp00.fits',
    'fastspec-loa-main-dark-nside1-hp01.fits',
    'fastspec-loa-main-dark-nside1-hp02.fits',
    'fastspec-loa-main-dark-nside1-hp03.fits',
    'fastspec-loa-main-dark-nside1-hp04.fits',
    'fastspec-loa-main-dark-nside1-hp05.fits',
    'fastspec-loa-main-dark-nside1-hp06.fits',
    'fastspec-loa-main-dark-nside1-hp07.fits',
    'fastspec-loa-main-dark-nside1-hp08.fits',
    'fastspec-loa-main-dark-nside1-hp09.fits',
    'fastspec-loa-main-dark-nside1-hp10.fits',
    'fastspec-loa-main-dark-nside1-hp11.fits',
    'fastspec-loa-special-backup.fits',
    'fastspec-loa-special-bright.fits',
    'fastspec-loa-special-dark.fits',
    'fastspec-loa-sv1-backup.fits',
    'fastspec-loa-sv1-bright.fits',
    'fastspec-loa-sv1-dark.fits',
    'fastspec-loa-sv1-other.fits',
    'fastspec-loa-sv2-backup.fits',
    'fastspec-loa-sv2-bright.fits',
    'fastspec-loa-sv2-dark.fits',
    'fastspec-loa-sv3-backup.fits',
    'fastspec-loa-sv3-bright.fits',
    'fastspec-loa-sv3-dark.fits' ]

specphot_catalogs = fastspec_catalogs
metadata_catalogs = fastspec_catalogs


## Define the columns that we want to extract from the fastspec catalog.  This can be edited as one needs. 
## The quantities available can be found in the fastspecfit data model: 
## https://fastspecfit.readthedocs.io/en/latest/fastspec.html

# cols_selection_metadata = [
#     'TARGETID', 'SURVEY', 'PROGRAM', 'DESI_TARGET', 'BGS_TARGET', 'SV1_BGS_TARGET', 'SV2_BGS_TARGET', 'SV3_BGS_TARGET', 'HEALPIX', 
#     'RA', 'DEC', 'Z', 'ZWARN', 'DELTACHI2', 'SPECTYPE', 'EBV', 
#     'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4', 
#     'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FLUX_IVAR_W3', 'FLUX_IVAR_W4']

# cols_selection_specphot = [
#     'RCHI2', 'RCHI2_LINE', 'RCHI2_CONT', 
#     'VDISP', 'VDISP_IVAR', 'TAUV', 'TAUV_IVAR', 'AGE', 'AGE_IVAR', 'ZZSUN', 'ZZSUN_IVAR', 
#     'DN4000', 'DN4000_OBS', 'DN4000_IVAR', 'DN4000_MODEL', 'DN4000_MODEL_IVAR', 
#     'LOGMSTAR', 'LOGMSTAR_IVAR', 'SFR', 'SFR_IVAR']

# cols_selection_fastspec = [
#     'APERCORR', 'APERCORR_R',
#     'INIT_BALMER_BROAD', 'INIT_SIGMA_NARROW', 'INIT_SIGMA_BALMER', 
#     'HALPHA_AMP', 'HALPHA_AMP_IVAR', 'HALPHA_FLUX', 'HALPHA_FLUX_IVAR', 'HALPHA_EW', 'HALPHA_SIGMA', 'HALPHA_SIGMA_IVAR',
#     'HALPHA_BROAD_AMP', 'HALPHA_BROAD_AMP_IVAR', 'HALPHA_BROAD_FLUX', 'HALPHA_BROAD_FLUX_IVAR', 'HALPHA_BROAD_SIGMA', 'HALPHA_BROAD_SIGMA_IVAR',
#     'HBETA_AMP', 'HBETA_AMP_IVAR', 'HBETA_FLUX', 'HBETA_FLUX_IVAR', 'HBETA_EW', 'HBETA_SIGMA', 'HBETA_SIGMA_IVAR',
#     'HBETA_BROAD_AMP', 'HBETA_BROAD_AMP_IVAR', 'HBETA_BROAD_FLUX', 'HBETA_BROAD_FLUX_IVAR', 'HBETA_BROAD_SIGMA', 'HBETA_BROAD_SIGMA_IVAR',
#     'OIII_5007_AMP', 'OIII_5007_AMP_IVAR', 'OIII_5007_FLUX', 'OIII_5007_FLUX_IVAR', 'OIII_5007_SIGMA', 'OIII_5007_SIGMA_IVAR', 
#     'OII_3726_AMP', 'OII_3726_AMP_IVAR', 'OII_3726_FLUX', 'OII_3726_FLUX_IVAR', 'OII_3726_SIGMA', 'OII_3726_SIGMA_IVAR', 
#     'OII_3729_AMP', 'OII_3729_AMP_IVAR', 'OII_3729_FLUX', 'OII_3729_FLUX_IVAR', 'OII_3729_SIGMA', 'OII_3729_SIGMA_IVAR', 
#     'NII_6548_AMP', 'NII_6548_AMP_IVAR', 'NII_6548_FLUX', 'NII_6548_FLUX_IVAR', 'NII_6548_SIGMA', 'NII_6548_SIGMA_IVAR', 
#     'NII_6584_AMP', 'NII_6584_AMP_IVAR', 'NII_6584_FLUX', 'NII_6584_FLUX_IVAR', 'NII_6584_SIGMA', 'NII_6584_SIGMA_IVAR', 
#     'OIII_5007_CONT', 'OIII_5007_CONT_IVAR', 'HALPHA_CONT', 'HALPHA_CONT_IVAR']

#############
# From manasvee

cols_selection_metadata = [
    'TARGETID', 'SURVEY', 'PROGRAM',
    'RA', 'DEC', # deg
    'BGS_TARGET', #'SV1_BGS_TARGET', 'SV2_BGS_TARGET', 'SV3_BGS_TARGET',
    'Z', 'ZWARN', 'DELTACHI2', 'SPECTYPE', # Redshift based on Redrock or QuasarNet (for QSO targets only).
    'FLUX_G', 'FLUX_R', 'FLUX_Z', # nmgy; Total g, r and z-band flux corrected for Galactic extinction.
    'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', # 1 / nmgy2
    ]

cols_selection_specphot = [
    'LOGMSTAR', 'LOGMSTAR_IVAR', # Msun; Logarithmic stellar mass (h=1.0, Chabrier+2003 initial mass function).
    'ABSMAG01_SDSS_U', 'ABSMAG01_SDSS_G', 'ABSMAG01_SDSS_R', 'ABSMAG01_SDSS_I', 'ABSMAG01_SDSS_Z', # Absolute magnitude in SDSS u-band band-shifted to z=0.1 assuming h=1.0.
    'ABSMAG01_IVAR_SDSS_U', 'ABSMAG01_IVAR_SDSS_G', 'ABSMAG01_IVAR_SDSS_R', 'ABSMAG01_IVAR_SDSS_I', 'ABSMAG01_IVAR_SDSS_Z',
    ] 

cols_selection_fastspec = [
    'APERCORR', 'APERCORR_R',
    'INIT_BALMER_BROAD', 'INIT_SIGMA_NARROW', 'INIT_SIGMA_BALMER', 
    'HALPHA_AMP', 'HALPHA_AMP_IVAR', 'HALPHA_FLUX', 'HALPHA_FLUX_IVAR', 'HALPHA_EW', 'HALPHA_SIGMA', 'HALPHA_SIGMA_IVAR',
    'HALPHA_BROAD_AMP', 'HALPHA_BROAD_AMP_IVAR', 'HALPHA_BROAD_FLUX', 'HALPHA_BROAD_FLUX_IVAR', 'HALPHA_BROAD_SIGMA', 'HALPHA_BROAD_SIGMA_IVAR',
    'HBETA_AMP', 'HBETA_AMP_IVAR', 'HBETA_FLUX', 'HBETA_FLUX_IVAR', 'HBETA_EW', 'HBETA_SIGMA', 'HBETA_SIGMA_IVAR',
    'HBETA_BROAD_AMP', 'HBETA_BROAD_AMP_IVAR', 'HBETA_BROAD_FLUX', 'HBETA_BROAD_FLUX_IVAR', 'HBETA_BROAD_SIGMA', 'HBETA_BROAD_SIGMA_IVAR',
    'OIII_5007_AMP', 'OIII_5007_AMP_IVAR', 'OIII_5007_FLUX', 'OIII_5007_FLUX_IVAR', 'OIII_5007_SIGMA', 'OIII_5007_SIGMA_IVAR', 
    'OII_3726_AMP', 'OII_3726_AMP_IVAR', 'OII_3726_FLUX', 'OII_3726_FLUX_IVAR', 'OII_3726_SIGMA', 'OII_3726_SIGMA_IVAR', 
    'OII_3729_AMP', 'OII_3729_AMP_IVAR', 'OII_3729_FLUX', 'OII_3729_FLUX_IVAR', 'OII_3729_SIGMA', 'OII_3729_SIGMA_IVAR', 
    'NII_6548_AMP', 'NII_6548_AMP_IVAR', 'NII_6548_FLUX', 'NII_6548_FLUX_IVAR', 'NII_6548_SIGMA', 'NII_6548_SIGMA_IVAR', 
    'NII_6584_AMP', 'NII_6584_AMP_IVAR', 'NII_6584_FLUX', 'NII_6584_FLUX_IVAR', 'NII_6584_SIGMA', 'NII_6584_SIGMA_IVAR', 
    'OIII_5007_CONT', 'OIII_5007_CONT_IVAR', 'HALPHA_CONT', 'HALPHA_CONT_IVAR']

#############

## Define the function which determines which galaxies to extract from the entire Loa sample. 
## In the example here, I select all galaxies which have SPECTYPE=GALAXY and z<0.05

def get_selected_rows_metadata(cat):
    '''Return boolean array of rows that meet criteria.'''
    no_row = len(cat)
    rows = (
        (cat['SPECTYPE'] == 'GALAXY') &    # Flag for galaxies
        (cat['Z'] >= 0.01) &                # Lower limit on redshift   
        (cat['Z'] <= 0.06)                  # Upper limit on redshift
        # (cat['DELTACHI2'] >= 25) &        # Lower limit on delta chi2: increases with purity of correctly identified redshifts
        # (cat['ZWARN'] == 0) &             # Flag for redshift/spectype fitting problems (ZWARN=0 is good)
        # (cat['BGS_TARGET'] != 0.)
    )

    return rows, no_row

## Now we go through all the fastspecfit sub-catalogs and identify the objects that match our selection criteria: 
print('Looking for galaxies matching the criteria in the  fastspecfit files')

all_rows = []
counter = 0
for i in range(len(fastspec_catalogs)):
    print('Processing ' + fastspec_catalogs[i])
    cat_meta = fitsio.read(fastspec_path+metadata_catalogs[i], 'METADATA', columns=cols_selection_metadata)
    rows, no_row = get_selected_rows_metadata(cat_meta)
    rows = np.arange(len(rows))[rows] # change to row numbers instead of boolean
    all_rows.append(rows)
    counter += no_row
print(f'Found {counter:,} galaxies matching the criteria in the fastspecfit catalogs')
print(len(all_rows), 'sub-catalogs processed')

## We can now loop through all the sub-catlogs, and extract the galaxies we want 

specphot_tables = []
fast_tables = []
meta_tables = []

for i in range(len(fastspec_catalogs)):
    cat = Table(fitsio.read(fastspec_path+specphot_catalogs[i], 'SPECPHOT', rows=all_rows[i]))
    specphot_tables.append(cat)

    cat_meta = Table(fitsio.read(fastspec_path+metadata_catalogs[i], 'METADATA', rows=all_rows[i]))
    meta_tables.append(cat_meta)

    cat_fast = Table(fitsio.read(fastspec_path+fastspec_catalogs[i], 'FASTSPEC', rows=all_rows[i]))
    fast_tables.append(cat_fast)

fast_combined = vstack(fast_tables)
specphot_combined = vstack(specphot_tables)
meta_combined = vstack(meta_tables)

## Apply the filter to keep only the columns we identified above

meta_combined.keep_columns(cols_selection_metadata)
specphot_combined.keep_columns(cols_selection_specphot)
fast_combined.keep_columns(cols_selection_fastspec)


## Combine all the columns extracted from the three different extensions into one catalog, and save as a fits table 

temp = hstack([meta_combined, specphot_combined])
merged_cat = hstack([temp, fast_combined])
merged_cat.write('./loa-combined-lowz.fits', overwrite=True)

cat = merged_cat

print('Finished with the base extraction')

## This could be the end of it, but I want to add extra columns to my catalog from two different additional sources 
## (1) the full redshift catalogs (to get the flags that are important to remove duplicates) 
## (2) the imaging catalogs (to get quantities related to the shape and size of the galaxies) 

## Let's start with the redshift catalog 
print('Matching with the resdhift catalog')

zcatfile = '/global/cfs/cdirs/desi/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits'

## Read only the redshift and spectype columns to be able to use my get_selected_rows_metadata function to determine 
## the row numbers of my objects of interest. (this should obviously be modified as need be, for whatever selection criteria are desired) 

zcat = fitsio.read(zcatfile, columns=['Z', 'SPECTYPE'])
rows, _ = get_selected_rows_metadata(zcat)
rows = np.arange(len(rows))[rows]  # transform into indices rather than boolean array

## Read in those rows, transform the result into an Astropy table, select which columns to keep, and
## merge those columns with the fastspecfit catalog 

zlowz = fitsio.read(zcatfile, rows=rows)
ztab = Table(zlowz)
ztab.keep_columns(['TARGETID', 'SURVEY', 'PROGRAM', 'DESINAME', 'MORPHTYPE', 'MAIN_PRIMARY', 'SV_PRIMARY', 'ZCAT_NSPEC', 'ZCAT_PRIMARY'])

cat_wflags = join(cat, ztab, join_type='left', keys=['TARGETID', 'SURVEY', 'PROGRAM'])


## Remove the duplicates by only keeping those with ZCAT_PRIMARY==True
mask = (cat_wflags['ZCAT_PRIMARY']  == True)
cat = cat_wflags[mask]

cat.write('loa-combined-lowz-zflags.fits', overwrite=True)

## Let's now add parameters from the Legacy Survey DR9 photometry catalog 
print('Matching with the Legacy Survey catalog')

specprod = 'loa'
vacdir = '/global/cfs/cdirs/desi/vac/dr2/lsdr9-photometry/loa/v1.0/observed-targets'

def read_tractorphot(cat, specprod='loa', verbose=False):
    """Gather Tractor photometry for an input catalog. Note that this function 
    hasn't been tested very extensively and will likely fail if there are duplicate
    TARGETIDs in the input redshift catalog.
    Note that some objects are missing the tractor photometry, and are therefore not included in the output 
    We should therefore not expect to have a perfect line-to-line match between the input and output catalogs of this function! 
    """
    from glob import glob
    from desitarget import geomask
    from desimodel.footprint import radec2pix

    tractorphotfiles = glob(os.path.join(vacdir, 'tractorphot', f'tractorphot-nside4-hp???-{specprod}.fits'))

    print(tractorphotfiles[0])
    hdr = fitsio.read_header(tractorphotfiles[0], 'TRACTORPHOT')
    tractorphot_nside = hdr['FILENSID']

    pixels = radec2pix(tractorphot_nside, cat['RA'], cat['DEC'])
    phot = []
    for pixel in sorted(set(pixels)):
        J = pixel == pixels
        photfile = os.path.join(vacdir, 'tractorphot', f'tractorphot-nside4-hp{pixel:03d}-{specprod}.fits')
        if os.path.isfile(photfile):
            targetids = fitsio.read(photfile, columns='TARGETID')
            K = np.where(np.isin(targetids, cat['TARGETID'][J]))[0]

            if verbose:
                print(f'Gathering Tractor photometry for {len(K):,d} object(s) from {photfile}')

            _phot = fitsio.read(photfile, rows=K)
            phot.append(Table(_phot))

    if len(phot) > 0:
        phot = vstack(phot)
    else:
        phot = Table()

    return phot


tractor = read_tractorphot(cat, specprod=specprod, verbose=True)
tractor

tractor.keep_columns(['TARGETID', 'TYPE', 'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'SHAPE_R', 'SHAPE_E1', 'SHAPE_E2', 'SHAPE_R_IVAR', 'SHAPE_E1_IVAR', 'SHAPE_E2_IVAR', 'SERSIC', 'SERSIC_IVAR'])

cat_withphot = join(cat, tractor, join_type='left', keys=['TARGETID'])
cat_withphot.write('./loa-combined-lowz-fastspec-phot.fits', overwrite=True)
print(len(cat_withphot), 'objects in the final catalog')
print('Done')
