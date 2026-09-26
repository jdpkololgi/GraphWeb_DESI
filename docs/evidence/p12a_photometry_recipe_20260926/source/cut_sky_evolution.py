import numpy as np
import matplotlib.pyplot as plt
import h5py

from hodpy.catalogue import Catalogue
from hodpy.luminosity_function import LuminosityFunctionTargetBGS
from hodpy.colour import ColourNew
from hodpy.hod_bgs import HOD_BGS_Simple
from hodpy import lookup


def cut_sky(position, velocity, magnitude, is_cen, cosmology, Lbox, zsnap, kcorr_r, kcorr_g,
            hod, replication=(0,0,0), zcut=None, mag_cut=None, cosmology_orig=None):
    """
    Creates a cut sky mock by converting the cartesian coordiantes of a cubic box mock to ra, dec, z.
    Magnitudes and colours are evolved with redshift
    Args:
        position:  array of comoving position vectors (Mpc/h), in the range -Lbox/2 < pos < Lbox/2
        velocity:  array of proper velocity vectors (km/s)
        magnitude: array of absolute magnitude
        is_cen:    boolean array indicating if galaxy is central (True) or satellite (False)
        cosmology: instance of astropy.cosmology class
        Lbox:      comoving box length of simulation (Mpc/h)
        zsnap:     redshift of simulation snapshot
        kcorr_r:   GAMA_KCorrection object with r-band k-correction
        kcorr_g:   GAMA_KCorrection object with g-band k-correction
        replication: tuple indicating which periodic replication to use. Default value is (0,0,0) 
                         (ie no replications).
        zcut:    If provided, will only return galaxies with z<=zcut. By default will return
                         all galaxies.
        mag_cut: If provided, will only return galaxies with apparent magnitude < mag_cut. 
                         By default will return all galaxies.
        cosmology_orig: instance of astropy.cosmology class. The original simulation cosmology.
                         If provided, magnitudes will be scaled by cosmology
    Returns:
        ra:   array of ra (deg)
        dec:  array of dec (deg)
        zcos: array of cosmological redshift, which does not include the effect of peculiar velocities
        zobs: array of observed redshift, which includes peculiar velocities.
        magnitude_new: array of new absolute magnitude, rescaled to match target luminosity 
                          function at each redshift
        app_mag: array of apparent magnitudes (calculated from rescaled magnitudes and colours)
        colour_new: array of g-r colours, which are re-assigned to add evolution
        colour_obs: array of observer-frame g-r colours
        index: array of indices. Used to match galaxies between the input and output arrays of 
                  this function
    """
    
    index = np.arange(position.shape[0])
    
    cat = Catalogue(cosmology)
    
    position_rep = position.copy()
    if replication==(0,0,0):
        print("No periodic replications")
    else:
        print("Applying periodic replications")
        for i in range(3):
            print("%.1f < %s < %.1f"%((-1+2*replication[i])*Lbox/2., chr(120+i), (1+2*replication[i])*Lbox/2.))
            position_rep[:,i] += Lbox*replication[i]
    
    ra, dec, zcos = cat.pos3d_to_equitorial(position_rep)
    vlos = cat.vel_to_vlos(position_rep, velocity)
    zobs = cat.vel_to_zobs(zcos, vlos)
    
    if not zcut is None:
        print("Applying redshift cut z < %.2f"%zcut)
        keep = zobs <= zcut
        ra, dec, zcos, zobs, magnitude, is_cen, index = \
                ra[keep], dec[keep], zcos[keep], zobs[keep], magnitude[keep], is_cen[keep], index[keep]
                  
    print("Rescaling magnitudes")
    lf = LuminosityFunctionTargetBGS(target_lf_file=lookup.target_lf, 
                                     sdss_lf_file=lookup.sdss_lf_tabulated, 
                                     lf_param_file=lookup.gama_lf_fits, 
                    hod_bgs_simple=HOD_BGS_Simple(lookup.bgs_hod_parameters))
    
    # first rescale to get target LF exactly
    #magnitude_new = lf.rescale_magnitude_to_target_box(magnitude, zsnap, volume,
    #                                cosmo_orig=cosmology_orig, cosmo_new=cosmology)
    # then rescale to get evolving target LF
    magnitude_new = lf.rescale_magnitude(magnitude, np.ones(len(zobs))*zsnap, zobs,
                                        cosmo_orig=cosmology_orig, cosmo_new=cosmology)
    
    print("Assigning colours")
    is_sat = np.invert(is_cen)
    colour_new = np.zeros(len(magnitude_new))
    
    #col = Colour()
    col = ColourNew(hod=hod)
    
    # randomly assign colours to centrals and satellites
    colour_new[is_cen] = col.get_central_colour(magnitude_new[is_cen], zobs[is_cen])
    colour_new[is_sat] = col.get_satellite_colour(magnitude_new[is_sat], zobs[is_sat])
    
    
    # get apparent magnitude
    app_mag = kcorr_r.apparent_magnitude(magnitude_new, zobs, colour_new)
    
    # observer frame colours
    colour_obs = colour_new + kcorr_g.k(zobs, colour_new) - kcorr_r.k(zobs, colour_new)
    
    if not mag_cut is None:
        print("Applying magnitude cut r < %.2f"%mag_cut)
        keep = app_mag <= mag_cut
        ra, dec, zcos, zobs, magnitude_new, app_mag, colour_new, colour_obs, index = \
                    ra[keep], dec[keep], zcos[keep], zobs[keep], magnitude_new[keep], \
                    app_mag[keep], colour_new[keep], colour_obs[keep], index[keep]
        
    
    return ra, dec, zcos, zobs, magnitude_new, app_mag, colour_new, colour_obs, index

