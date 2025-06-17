import astropy.constants as c
import astropy.units as u
from astropy.table import Table, join, vstack, unique, hstack, unique
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo
from astropy.io import fits
import fitsio
import os
import numpy as np

import matplotlib.pyplot as plt

class GalaxyCatalog:
    def __init__(self, PATH, LOGMSTAR=9.):
        self.zcat = Table(fitsio.read(PATH))
        # self.zcat = Table(fitsio.read(f'{specprod_dir}/zcatalog/v1/zpix-main-bright.fits', "ZCATALOG", columns=columns))
        self.zcat = self.zcat[self.zcat['ZWARN'] == 0] # Filter for successful redshifts
        self.zcat = self.zcat[self.zcat['DELTACHI2'] >= 25.] # Filter for good quality redshifts
        self.zcat = self.zcat[self.zcat['LOGMSTAR'] >= LOGMSTAR] # Filter for galaxies with M_STAR > LOGMSTAR
        self.zcat = self.zcat[self.zcat['BGS_TARGET'] != 0.] # Filter for BGS targets

        print(f"Raw number of zcat entries: {len(self.zcat):,}")

    def select_gal_classes(self, gal_class='BGS', M_STAR=9):
        '''
        Select galaxies based on their class.
        '''
        if gal_class not in ['BGS', 'ELG', 'LRG', 'QSO']:
            raise ValueError("gal_class must be one of ['BGS', 'ELG', 'LRG', 'QSO']")

        desi_tgt = self.zcat['DESI_TARGET']

        is_bgs  = (desi_tgt & 2**60 != 0)
        is_lrg  = (desi_tgt & 2**0 != 0)
        is_elg  = (desi_tgt & 2**1 != 0)
        is_qso  = (desi_tgt & 2**2 != 0)
        is_mws  = (desi_tgt & 2**61 != 0)
        is_scnd = (desi_tgt & 2**62 != 0)

        # Number of sources of each target type
        n_bgs = np.count_nonzero(is_bgs)
        n_lrg = np.count_nonzero(is_lrg)
        n_elg = np.count_nonzero(is_elg)
        n_qso = np.count_nonzero(is_qso)
        n_mws = np.count_nonzero(is_mws)
        n_scnd = np.count_nonzero(is_scnd)
        
        # Let us look at the numbers visually - 

        plt.figure(figsize = (8, 4))

        targets = ['BGS', 'LRG', 'ELG', 'QSO', 'MWS', 'SCND']
        numbers = [n_bgs, n_lrg, n_elg, n_qso, n_mws, n_scnd]

        plt.bar(targets, numbers, color = 'purple', alpha = 0.5)
        plt.ylabel('Number of primary spectra')
        plt.yscale('log')
        plt.tight_layout()

        if gal_class == 'BGS':
            self.zcat = self.zcat[is_bgs]
        elif gal_class == 'LRG':
            self.zcat = self.zcat[is_lrg]
        elif gal_class == 'ELG':
            self.zcat = self.zcat[is_elg]
        elif gal_class == 'QSO':
            self.zcat = self.zcat[is_qso]
        else:
            raise ValueError("gal_class must be one of ['BGS', 'ELG', 'LRG', 'QSO']")
        
        # Filter for galaxies with M_STAR > 9
        self.zcat = self.zcat[self.zcat['LOGMSTAR'] > M_STAR]

    def sky_coord(self, plot=True, globeplot=True):
        '''
        Return the SkyCoord object for the galaxies in the catalogue.
        '''

        if plot:
            if globeplot==False:
                plt.figure(figsize=(18, 12))
                plt.scatter(self.zcat['RA'], self.zcat['DEC'], c= self.zcat['Z'], cmap='plasma', s=0.1, alpha=0.5, marker='.')
                plt.xlim(0, 360)
                plt.xlabel('RA (degrees)')
                plt.ylabel('DEC (degrees)')
                plt.title(f'Sky Coordinates ({len(self.zcat):,} galaxies)')
                plt.colorbar(label='Redshift (z)')
                plt.grid()
                plt.show()
            else:
                # project coordinates onto a sky map
                fig = plt.figure(figsize=(18, 6))
                ax = fig.add_subplot(111, projection='mollweide')
                ax.scatter(np.radians(self.zcat['RA'])-np.pi, np.radians(self.zcat['DEC']), c=self.zcat['Z'], cmap='inferno_r', s=0.1, alpha=0.5, marker='.')
                # ax.set_xlim(0, 360)
                ax.set_xlabel('RA (rad)')
                ax.set_ylabel('DEC (rad)')
                ax.set_title(f'Sky Coordinates ({len(self.zcat):,} galaxies)')
                # plt.colorbar(label='Redshift (z)', ax=ax)
                plt.grid()
                plt.show()


        self.sky_coord = SkyCoord(ra=self.zcat['RA'], dec=self.zcat['DEC'], unit=(u.deg, u.deg))

    def cartesian_coord(self):
        '''
        Return the Cartesian coordinates of the galaxies in the catalogue assuming a Planck 2018 cosmology.
        '''
        self.comoving_distance = cosmo.comoving_distance(self.zcat['Z']).to(u.Mpc)
        self.cartesian_coord = SkyCoord(ra=self.zcat['RA'], dec=self.zcat['DEC'], unit=(u.deg, u.deg), distance=self.comoving_distance)

        self.X = self.cartesian_coord.cartesian.x.to(u.Mpc).value
        self.Y = self.cartesian_coord.cartesian.y.to(u.Mpc).value
        self.Z = self.cartesian_coord.cartesian.z.to(u.Mpc).value

        plt.style.use('dark_background')
        plt.rcParams.update({'font.size': 14})
        
        fig = plt.figure(figsize=(18, 12))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.X, self.Y, self.Z, c=self.zcat['Z'], cmap='inferno_r', s=0.1, alpha=0.5, marker='.')
        ax.set_xlabel('X (Mpc)')
        ax.set_ylabel('Y (Mpc)')
        ax.set_zlabel('Z (Mpc)')
        ax.grid(False)
        # zoom and orient the plot

        ax.view_init(elev=30, azim=60)
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.set_title(f'DESI_BGS_cartesian_10_9_dark ({len(self.zcat):,} galaxies)')
        plt.colorbar(ax.collections[0], label='Redshift (z)')
        fig.savefig('DESI_BGS_cartesian_10_9_dark.pdf')
        plt.show()

        
if __name__ == "__main__":
    # PATH = '/global/homes/d/dkololgi/GraphWeb_DESI/test.fits'
    PATH='/global/homes/d/dkololgi/GraphWeb_DESI/loa-combined-lowz.fits'
    data = GalaxyCatalog(PATH)
    # data.select_gal_classes('BGS') # Example usage to select BGS galaxies
    data.sky_coord(plot=True, globeplot=False) # Example usage to plot sky coordinates
    data.cartesian_coord() # Example usage to compute Cartesian coordinates