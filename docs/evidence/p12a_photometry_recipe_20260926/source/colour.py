#! /usr/bin/env python
import numpy as np
from scipy.special import erfc
from scipy.interpolate import RegularGridInterpolator, interp1d
from hodpy import lookup


class ColourNew(object):
    
    def __init__(self, colour_fits=lookup.colour_fits, hod=None):
        """
        Class containing methods for randomly assigning galaxies a g-r colour
        using updated colour distributions
        """
        
        self.hod = hod
        
        self.redshift_bin, self.redshift_median, \
                self.functions, self.parameters = self.read_fits(colour_fits)
        
        
        self.__mag_bins = np.arange(-23, -14, 0.01)
        
        self.__blue_mean_interpolator = self.__get_interpolator(0)
        self.__blue_rms_interpolator = self.__get_interpolator(1)
        self.__red_mean_interpolator = self.__get_interpolator(2)
        self.__red_rms_interpolator = self.__get_interpolator(3)
        self.__fraction_blue_interpolator = self.__get_interpolator(4)
        
        self.__central_fraction_interpolator = self.__initialize_central_fraction_interpolator(z=0.2)
            
            
    def __initialize_central_fraction_interpolator(self, z=0.2):
        
        magnitudes = np.arange(-23,-10,0.1)
        fcen = np.zeros(len(magnitudes))

        for i in range(len(magnitudes)):
            magnitude = np.array([magnitudes[i],])
            redshift = np.array([z,])
            logMmin = np.log10(self.hod.Mmin(magnitude))
            logM1 = np.log10(self.hod.M1(magnitude))
            logM0 = np.log10(self.hod.M0(magnitude))
            sigmalogM = self.hod.sigma_logM(magnitude)
            alpha = self.hod.alpha(magnitude)

            n_all = self.hod.get_n_HOD(magnitude, redshift, logMmin, logM1, logM0, sigmalogM, alpha,
                            Mmin=10, Mmax=16, galaxies="all")

            n_cen = self.hod.get_n_HOD(magnitude, redshift, logMmin, logM1, logM0, sigmalogM, alpha,
                            Mmin=10, Mmax=16, galaxies="cen")

            fcen[i] = n_cen/n_all
            
        magnitudes2 = np.arange(-28,10,0.1)
        fcen2 = np.zeros(len(magnitudes2))
        fcen2[50:180] = fcen
        fcen2[:50] = fcen[0]
        fcen2[180:] = fcen[-1]
    
        return interp1d(magnitudes2, fcen2, kind='cubic')

        
        
    def read_fits(self, colour_fits):
        fits = np.load(colour_fits)
        Nbins = fits.shape[0] # number of redshift bins
        redshift_bin    = np.zeros(Nbins) # bin centres
        redshift_median = np.zeros(Nbins) # median redshift in bin
        functions       = np.zeros((5,Nbins),dtype="i")
        parameters      = [None]*Nbins

        for i in range(Nbins):
            redshift_bin[i] = fits[i,0,0]
            redshift_median[i] = fits[i,0,1]
            functions[:,i] = fits[i,1:,0]
            parameters[i] = fits[i,1:,1:]
        
        return redshift_bin, redshift_median, functions, parameters
            
        
    def broken(self, x, a, b, c, d):
        """
        Broken linear function with smooth transition
        """
        trans=20
        y1 = a*x + b
        y2 = c*x + d
        return np.log10(10**((y1)*trans) + 10**((y2)*trans)) / trans


    def broken_reverse(self, x, a, b, c, d):
        """
        Broken linear function with smooth transition
        """
        trans=20
        y1 = a*x + b
        y2 = c*x + d
        return 1-np.log10(10**((1-y1)*trans) + 10**((1-y2)*trans)) / trans
    
    
    def __get_interpolator(self, param_idx):
        
        redshifts = self.redshift_bin.copy()
        params = self.parameters.copy()
            
        array = np.zeros((len(self.__mag_bins), len(redshifts)))
        
        for i in range(len(redshifts)):
            if self.functions[param_idx][i] == 0:
                array[:,i] = self.broken(self.__mag_bins, *params[i][param_idx])
            else:
                array[:,i] = self.broken_reverse(self.__mag_bins, *params[i][param_idx])
                
        array = np.clip(array, -10, 10)
        
        func = RegularGridInterpolator((self.__mag_bins, redshifts), array,
                                       method='linear', bounds_error=False, fill_value=None)
        return func
        
        
    def blue_mean(self, magnitude, redshift):
        return self.__blue_mean_interpolator((magnitude, redshift))
    
    def blue_rms(self, magnitude, redshift):
        return np.clip(self.__blue_rms_interpolator((magnitude, redshift)), 0.02, 10)
    
    def red_mean(self, magnitude, redshift):
        return self.__red_mean_interpolator((magnitude, redshift))
    
    def red_rms(self, magnitude, redshift):
        return np.clip(self.__red_rms_interpolator((magnitude, redshift)), 0.02, 10)
        
    def fraction_blue(self, magnitude, redshift):
        frac_blue = np.clip(self.__fraction_blue_interpolator((magnitude, redshift)), 0, 1)
        
        # if at bright end blue_mean > red_mean, set all galaxies as being red
        b_m = self.blue_mean(magnitude, redshift)
        r_m = self.red_mean(magnitude, redshift)
        idx = np.logical_and(b_m > r_m, magnitude<-20)
        frac_blue[idx] = 0
        
        return frac_blue
    
    
    def fraction_central(self, magnitude, redshift):
        """
        Fraction of central galaxies as a function of magnitude and redshift

        Args:
            magnitude: array of absolute r-band magnitudes (with h=1)
            redshift:  array of redshifts
        Returns:
            array of fraction of central galaxies
        """
        
        # for Abacus
        # nsat_ncen = 0.509 * (2 - erfc(0.4467*(magnitude+20.26)))
        # return 1 / (1 + nsat_ncen)
        
        return self.__central_fraction_interpolator(magnitude)
        
        
        
    def satellite_mean(self, magnitude, redshift):
        """
        Mean satellite colour as a function of magnitude and redshift

        Args:
            magnitude: array of absolute r-band magnitudes (with h=1)
            redshift:  array of redshifts
        Returns:
            array of g-r colours
        """
        
        #colour = 0.86 - 0.065 * (magnitude + 20)
        #ind = redshift > 0.1
        #colour[ind] -= 0.18 * (redshift[ind]-0.1) 

        blue_mean = self.blue_mean(magnitude, redshift)
        red_mean = self.red_mean(magnitude, redshift)
        
        colour = (0.8*red_mean + 0.2*blue_mean)
        
        return colour 
    
    
    def probability_red_satellite(self, magnitude, redshift):
        """
        Probability a satellite is red as a function of magnitude and redshift

        Args:
            magnitude: array of absolute r-band magnitudes (with h=1)
            redshift:  array of redshifts
        Returns:
            array of probabilities
        """
        
        sat_mean  = self.satellite_mean(magnitude, redshift)
        blue_mean = self.blue_mean(magnitude, redshift)
        red_mean  = self.red_mean(magnitude, redshift)
        
        #p_red = np.clip(np.absolute(sat_mean-blue_mean) / \
        #                np.absolute(red_mean-blue_mean), 0, 1)
        p_red = np.clip((sat_mean-blue_mean) / \
                        (red_mean-blue_mean), 0, 1)
        f_blue = self.fraction_blue(magnitude, redshift)
        
        idx = f_blue==0
        p_red[idx]=1
        idx = f_blue==1
        p_red[idx]=0
        
        f_cen = self.fraction_central(magnitude, redshift)

        return np.minimum(p_red, ((1-f_blue)/(1-f_cen)))
    
    
    def probability_red_central(self, magnitude, redshift):
        prob_red_sat  = self.probability_red_satellite(magnitude, redshift)
        prob_blue_sat = 1. - prob_red_sat

        frac_cent = self.fraction_central(magnitude, redshift)
        frac_blue = self.fraction_blue(magnitude, redshift)

        prob_blue = frac_blue/frac_cent - prob_blue_sat/frac_cent + prob_blue_sat
        
        return 1. - prob_blue
    
    
    
    
    def get_satellite_colour(self, magnitude, redshift):
        """
        Randomly assigns a satellite galaxy a g-r colour

        Args:
            magnitude: array of absolute r-band magnitudes (with h=1)
            redshift:  array of redshifts
        Returns:
            array of g-r colours
        """

        num_galaxies = len(magnitude)

        # probability the satellite should be drawn from the red sequence
        prob_red = self.probability_red_satellite(magnitude, redshift)

        # random number for each galaxy 0 <= u < 1
        u = np.random.rand(num_galaxies)

        # if u <= p_red, draw from red sequence, else draw from blue sequence
        is_red = u <= prob_red
        is_blue = np.invert(is_red)
    
        mean = np.zeros(num_galaxies, dtype="f")
        mean[is_red]  = self.red_mean(magnitude[is_red],   redshift[is_red])
        mean[is_blue] = self.blue_mean(magnitude[is_blue], redshift[is_blue])

        stdev = np.zeros(num_galaxies, dtype="f")
        stdev[is_red]  = self.red_rms(magnitude[is_red],   redshift[is_red])
        stdev[is_blue] = self.blue_rms(magnitude[is_blue], redshift[is_blue])

        # randomly select colour from Gaussian
        colour = np.random.normal(loc=0.0, scale=1.0, size=num_galaxies)
        colour = colour * stdev + mean

        return colour


    def get_central_colour(self, magnitude, redshift):
        """
        Randomly assigns a central galaxy a g-r colour

        Args:
            magnitude: array of absolute r-band magnitudes (with h=1)
            redshift:  array of redshifts
        Returns:
            array of g-r colours
        """
        num_galaxies = len(magnitude)

        # find probability the central should be drawn from the red sequence
        prob_red_sat  = self.probability_red_satellite(magnitude, redshift)
        prob_blue_sat = 1. - prob_red_sat

        frac_cent = self.fraction_central(magnitude, redshift)
        frac_blue = self.fraction_blue(magnitude, redshift)

        prob_blue = frac_blue/frac_cent - prob_blue_sat/frac_cent + \
                                                          prob_blue_sat
        prob_red = 1. - prob_blue

        # random number for each galaxy 0 <= u < 1
        u = np.random.rand(num_galaxies)

        # if u <= p_red, draw from red sequence, else draw from blue sequence
        is_red = u <= prob_red
        is_blue = np.invert(is_red)

        mean = np.zeros(num_galaxies, dtype="f")
        mean[is_red]  = self.red_mean(magnitude[is_red],   redshift[is_red])
        mean[is_blue] = self.blue_mean(magnitude[is_blue], redshift[is_blue])

        stdev = np.zeros(num_galaxies, dtype="f")
        stdev[is_red]  = self.red_rms(magnitude[is_red],   redshift[is_red])
        stdev[is_blue] = self.blue_rms(magnitude[is_blue], redshift[is_blue])

        # randomly select colour from gaussian
        colour = np.random.normal(loc=0.0, scale=1.0, size=num_galaxies)
        colour = colour * stdev + mean

        return colour
