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
        self.cartesian_coord() # Generate Cartesian coordinates for north and south on initialization
        print(f'Reordering zcat to match [north, south] Cartesian coordinates...')
        self.zcat_north = self.zcat[self.galactic_north_mask].copy()
        self.zcat_south = self.zcat[self.galactic_south_mask].copy()
        self.zcat = vstack([self.zcat_north, self.zcat_south])
        print(f"Number of zcat entries after filtering: {len(self.zcat):,}")
    
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

    def compute_sky_coord(self, plot=True, globeplot=True):
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

        self.comoving_distance = cosmo.comoving_distance(self.zcat['Z']).to(u.Mpc)

        self.sky_coord_icrs = SkyCoord(
            ra=self.zcat['RA'],
            dec=self.zcat['DEC'], 
            unit=(u.deg, u.deg), 
            distance=self.comoving_distance,
            frame='icrs'
        )
        # Convert to Galactic coordinates
        sky_gal = self.sky_coord_icrs.galactic
        self.galactic_north_mask = sky_gal.b.deg > 0.
        self.galactic_south_mask = ~self.galactic_north_mask

        # self.sky_coord_galactic = self.sky_coord.galactic
        # self.galactic_north_mask = self.sky_coord_galactic.b.deg > 0.
        # self.galactic_south_mask = self.sky_coord_galactic.b.deg < 0.
        # self.sky_coord_galactic_north = self.sky_coord_galactic[self.galactic_north_mask]
        # self.sky_coord_galactic_south = self.sky_coord_galactic[self.galactic_south_mask]
        
    def cartesian_coord(self, xyzplot='hemispheres'):
        '''
        Return the Cartesian coordinates of the galaxies in the catalogue assuming a Planck 2018 cosmology.
        '''
        self.compute_sky_coord(plot=False, globeplot=False)  # Ensure sky coordinates are computed first

        # get full cartesian coordinates
        cart = self.sky_coord_icrs.cartesian
        self.X = cart.x.to(u.Mpc).value
        self.Xn = cart[self.galactic_north_mask].x.to(u.Mpc).value
        self.Xs = cart[self.galactic_south_mask].x.to(u.Mpc).value

        self.Y = cart.y.to(u.Mpc).value
        self.Yn = cart[self.galactic_north_mask].y.to(u.Mpc).value
        self.Ys = cart[self.galactic_south_mask].y.to(u.Mpc).value

        self.Z = cart.z.to(u.Mpc).value
        self.Zn = cart[self.galactic_north_mask].z.to(u.Mpc).value
        self.Zs = cart[self.galactic_south_mask].z.to(u.Mpc).value

        # self.X = self.sky_coord.cartesian.x.to(u.Mpc).value
        # self.Xn = self.sky_coord_galactic_north.cartesian.x.to(u.Mpc).value
        # self.Xs = self.sky_coord_galactic_south.cartesian.x.to(u.Mpc).value
        # self.Y = self.sky_coord.cartesian.y.to(u.Mpc).value
        # self.Yn = self.sky_coord_galactic_north.cartesian.y.to(u.Mpc).value
        # self.Ys = self.sky_coord_galactic_south.cartesian.y.to(u.Mpc).value
        # self.Z = self.sky_coord.cartesian.z.to(u.Mpc).value
        # self.Zn = self.sky_coord_galactic_north.cartesian.z.to(u.Mpc).value
        # self.Zs = self.sky_coord_galactic_south.cartesian.z.to(u.Mpc).value

        if xyzplot == 'hemispheres':
            redshifts = np.concatenate((self.zcat['Z'][self.galactic_north_mask], self.zcat['Z'][self.galactic_south_mask]))
            X = np.concatenate((self.Xn, self.Xs))
            Y = np.concatenate((self.Yn, self.Ys))
            Z = np.concatenate((self.Zn, self.Zs))
        elif xyzplot == 'full':
            redshifts = self.zcat['Z']
            X = self.X
            Y = self.Y
            Z = self.Z
        else:
            raise ValueError("xyzplot must be one of ['hemispheres', 'full']")

        plt.style.use('seaborn-v0_8-darkgrid')
        plt.rcParams["text.usetex"] = False
        plt.rcParams.update({'font.size': 14})
        
        fig = plt.figure(figsize=(18, 12))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X, Y, Z, c=redshifts, cmap='inferno_r', s=0.1, alpha=0.5, marker='.')
        ax.set_xlim([-300, 300])
        ax.set_ylim([-300, 300])
        ax.set_zlim([-300, 300])
        ax.set_xlabel('X (Mpc)')
        ax.set_ylabel('Y (Mpc)')
        ax.set_zlabel('Z (Mpc)')
        ax.grid(False)
        # zoom and orient the plot

        ax.view_init(elev=30, azim=60)
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.set_title(f'DESI_BGS_cartesian_10_9_dark_{xyzplot} ({len(self.zcat):,} galaxies)')
        plt.colorbar(ax.collections[0], label='Redshift (z)')
        fig.savefig('DESI_BGS_cartesian_10_9_dark.pdf')
        plt.show()

        
if __name__ == "__main__":
    # PATH = '/global/homes/d/dkololgi/GraphWeb_DESI/test.fits'
    PATH='/global/homes/d/dkololgi/GraphWeb_DESI/loa-combined-lowz.fits'
    data = GalaxyCatalog(PATH)
    # data.select_gal_classes('BGS') # Example usage to select BGS galaxies
    # data.sky_coord(plot=True, globeplot=False) # Example usage to plot sky coordinates
    # data.cartesian_coord(xyzplot='hemispheres') # Example usage to compute Cartesian coordinates