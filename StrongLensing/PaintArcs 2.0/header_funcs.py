#!/usr/bin/env python
#-*- coding:utf-8 -*-
# ===================================================
# Authors:
# Carlos Brandt - 
# ===================================================


""" Module to deal with FITS (image) header entries """

import sys;
import math as m;
import astropy.io.fits as pyfits;
import numpy

try:
    import pywcs
except ImportError:
    import astropy.wcs as pywcs

def get_pixelscale(hdr):
    """ Read out header dimpix (pixel_scale) from CD_* (or CDELT*) keys

    Input:
     - hdr  <pyfits.Header> : Image header instance

    Output:
     - pixel_scale  float : Image sky-to-pixel scale ("/px)

    ---
    """

    try:
        pixel_scale = hdr['PIXSCALE'];
        return pixel_scale;
    except:
        pass;
        
    try:
        CD1_1 = float(hdr['CD1_1']);
        CD1_2 = float(hdr['CD1_2']);
        CD2_1 = float(hdr['CD2_1']);
        CD2_2 = float(hdr['CD2_2']);
    except:
        CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1;

    if ( CD1_1 == 1 and CD1_2 == 1 and CD2_1 == 1 and CD2_2 == 1 ):
        try:
            CD1_1 = float(hdr['CDELT1']);
            CD1_2 = 0.;
            CD2_1 = 0.;
            CD2_2 = float(hdr['CDELT2']);
        except:
            CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1;
            print("Unable to find the pixel size value");

    try:
        CUNIT1 = str(hdr['CUNIT1']);
        CUNIT2 = str(hdr['CUNIT2']);
    except:
        print >> sys.stderr, "Warning: Unable to find 'CUNIT[1|2]' parameters in header instance for image coordinates unit.";
        print >> sys.stderr, "         Degrees ('deg') is being used as default unit value.";
        CUNIT1 = 'deg';
        CUNIT2 = 'deg';

    conv_factor = 1.;
    CUNIT = CUNIT1
    # Convertion from degrees to arcsec units..
    #
    '''if ( CUNIT1 != CUNIT2 ):
        >> sys.stderr, print("Error: I do not have support for different per-axis pixel resolutions.");
	>> sys.stderr, print("Error: CUNIT1: %s , CUNIT2: %s") % (CUNIT1,CUNIT2)
	return None'''

    # Scale is given, using CD1_1, CD1_2, CD2_1, CD2_2 header keys:
    #
    pixel_scale = m.sqrt((CD1_1**2 + CD1_2**2 + CD2_1**2 + CD2_2**2)/2.) #* conv_factor;

    return pixel_scale;


def update_coordinates(hdr, x_ini, y_ini):
    """Update header information regarding world coordinate system"""


    LTV1 = None;
    LTV2 = None;

    NAXIS1 = int(hdr['NAXIS1']);
    NAXIS2 = int(hdr['NAXIS2']);

    try:
        CRPIX1 = float(hdr['CRPIX1']);
        CRPIX2 = float(hdr['CRPIX2']);
    except:
        CRPIX1 = CRPIX2 = 1.0;

    try:
        CD1_1 = float(hdr['CD1_1']);
        CD1_2 = float(hdr['CD1_2']);
        CD2_1 = float(hdr['CD2_1']);
        CD2_2 = float(hdr['CD2_2']);
    except:
        CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1.0;

    # Check the edges and update header keyworks for new image..
    #
    #	if ( x_ini >= 0  and  y_ini >= 0 ):
    #		LTV1 = -1*x_ini;
    #		LTV2 = -1*y_ini;
    #		CRPIX1 = CRPIX1 + LTV1;
    #		CRPIX2 = CRPIX2 + LTV2;
    #	if ( x_ini < 0  and  y_ini >= 0 ):
    #		LTV2 = -1*y_ini;
    #		CRPIX2 = CRPIX2 + LTV2;
    #	if ( y_ini < 0  and  x_ini >= 0 ):
    #		LTV1 = -1*x_ini;
    #		CRPIX1 = CRPIX1 + LTV1;
    LTV1 = -1*x_ini;
    LTV2 = -1*y_ini;
    CRPIX1 = CRPIX1 + LTV1;
    CRPIX2 = CRPIX2 + LTV2;


    # Add some keys to header..
    #
    WCSDIM = 2
    CDELT1  =  CD1_1;
    CDELT2  =  CD2_2;
    LTM1_1 = 1.0
    LTM2_2 = 1.0
    WAT0_001= 'system=image'
    WAT1_001= 'wtype=tan axtype=ra'
    WAT2_001= 'wtype=tan axtype=dec'


    # Header update..
    #
    if (LTV1 != None) :
        hdr.update(LTV1 = LTV1)
        hdr.update(CRPIX1 = CRPIX1)
    if (LTV2 != None) :
        hdr.update(LTV2 = LTV2)
        hdr.update(CRPIX2 = CRPIX2)
    hdr.update(WCSDIM = WCSDIM)
    hdr.update(CDELT1 = CDELT1)
    hdr.update(CDELT2 = CDELT2)
    hdr.update(LTM1_1 = LTM1_1)
    hdr.update(LTM2_2 = LTM2_2)
    hdr.update(WAT0_001 = WAT0_001)
    hdr.update(WAT1_001 = WAT1_001)
    hdr.update(WAT2_001 = WAT2_001)


    return (hdr);


def get_ra_dec_limits(image):

    """
    This function reads RA and DEC values correponding to the spatial limits of the sky image. 

    Input:
    - image <str>: name of image(tile) fits file.

    Output:
    - ra,dec <tuple>: two tuples with image sky region limits.
    """

    hdr = pyfits.getheader(image)
    
    wcs = pywcs.WCS(hdr) 

    naxis1 = hdr['NAXIS1']
    naxis2 = hdr['NAXIS2']
    ra_centre = hdr['CRVAL1']
    dec_centre = hdr['CRVAL2']
    pixcrd = numpy.array([[1,1],[naxis1,naxis2]])

    skycrd = wcs.wcs_pix2sky(pixcrd,1)

    ra_min = skycrd[0][0]
    dec_min = skycrd[0][1]
    ra_max = skycrd[1][0]
    dec_max = skycrd[1][1]


    return ((ra_min,ra_max),(dec_min,dec_max))
    
    
def get_header_parameter( image_file, *parargs ):
    """
    Read parameter value from image header. 

    Input:
    - image_file <str>: image FITS filename
    - parargs <str>: parameters

    Output:
    - parvalue <str>: header parameter value
    """

    _header = pyfits.getheader( image_file );

    # Get image parameters in header.
    #
    param_list = [];


    for _param in parargs :
        try:
            param_list.append( _header[_param] ); 

        except:
            #print( >> sys.stderr, "Warning: Unable to find this parameter in header instance for image coordinates unit.");
            param_list.append(None) 
       
    return (param_list);
    
