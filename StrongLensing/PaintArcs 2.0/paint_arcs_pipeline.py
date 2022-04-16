#!/usr/bin/python
# =====================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# Eduardo Valadão - eduardovaladao98@gmail.com
# =====================================================

""" This is the pipeline of the PaintArcs 2.0, it generates a FITS file with images of a simulated
source-galaxy lensed by a physical lens model chosen by the user. It also adds the effect of seeing
and Poisson noise to the lensed images in a background with characteristic noise of ground-based
telescopes. This is a work in progress, so my apologies to the user if any errors arise and if the
code is too slow. """
 
import math
import convolve 
import add_noise
import numpy as np 
import paint_arcs_v2 as pa
import astropy.io.fits as pyfits
from convert_arcsec_pix import convert_arcsec_pix

def run(params, seeing, background, pix_scale, suffix_name=''):
    ''' This function generates a FITS file of the chosen lens model with elliptical source
    based on the input parameter list called params, it also adds seeing and noise. '''

    # Creates the FITS file with the resulting lensed images, i.e. a header and an array of pixels: 
    hdu = pa.paint_arcs(params, pix_scale)
    arc_array = hdu.data
    arc_header = hdu.header
    pyfits.writeto("arc"+suffix_name+".fits", arc_array, arc_header, overwrite=True) 
    dim_x = arc_header["NAXIS1"]
    dim_y = arc_header["NAXIS2"]
    
    # Convolves the lensed images with a gaussian PSF to generate the seeing effect:
    seeing = convert_arcsec_pix(seeing, pix_scale)	
    sigma = seeing/(2.0*math.sqrt(2.0*math.log(2.0)))
    arc_header['APSF'] = (seeing, 'Arcs PSF (arcsec)')	
    arc_conv_array = convolve.gauss_convolution_fft(arc_array, 3.0, sigma)
    pyfits.writeto("arc_conv"+suffix_name+".fits", arc_conv_array, arc_header, overwrite=True)
        
    # Adds Poisson noise to the lensed images:
    arc_conv_noise_array = add_noise.add_poisson_noise(arc_conv_array)
    pyfits.writeto("arc_conv_noise"+suffix_name+".fits", arc_conv_noise_array, arc_header, overwrite=True)	
    
    # Creates a background image with Poisson noise:
    bkg_array = background*np.ones((dim_y,dim_x)) 
    bkg_noise_array = add_noise.add_poisson_noise(bkg_array)
    pyfits.writeto("background"+suffix_name+".fits", bkg_noise_array, overwrite=True)
    
    # Adds the convoluted lensed images with noise to the background image:
    arc_conv_noise_bkg = arc_conv_noise_array + bkg_noise_array
    arc_header['BKG'] = (background, 'Background Mean Value')	
    
    pyfits.writeto("arc_conv_noise_bkg"+suffix_name+".fits", arc_conv_noise_bkg, arc_header, overwrite=True)	


#==================================================================================================================
# params = [effective radius of the source (r0), x position of the source's center in arcseconds (s1),
# y position of the source's center in arcseconds (s2), angular orientation of the source in degrees (phie), 
# index of the Sersic profile (n), effective radius of the Sersic profile (r_e), ellipticity of the source (e), 
# magnitude of the lensed images (mag), magnitude of point zero (mag_zpt), ellipticity of the mass distribution 
# in the case of SIE and of the potential in the case of SIEP (f), external convergence (k), external shear (u), 
# value of the angular einstein radius in arcseconds (theta_e), x position of the image's center in arcseconds 
# (x1_0, I usually consider it 4), y position of the image's center in arcseconds (x2_0, I usually consider it 4), 
# value ks <= 1 for the convergence of the NFW lens, number associated with the model chosen by the user (model = 1 
# for ellipse, 2 for SIS, 3 for Point Lens, 4 for SIE and 5 for SIEP)] 
# 
# The quantities x1_0 and x2_0 are just a change in coordinates and they cannot be zero, that's why I usually 
# take x1_0 = x2_0 = 4. This is a problem in the part of pixelation of my code, so I apologize for this 
# inconvenience. The user could, for example, try params = [0.6, 4.75, 4.0, 90, 1, 0.15, 0.8, 20, 31.83, 
# 0, 0, 0, 4, 4, 4, 2] to verify if the code is working.
#==================================================================================================================

params = [0.6, 4.6, 4.6, 80, 1, 0.15, 0.8, 16, 31.83, 0, 0.2, 0.1, 4, 4, 4, 1.0, 6]
run(params, 0.5, 75, 0.04, "einstein")

