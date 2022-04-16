
""" This package contains functions convert from arcsec units to pixels and vice-versa. """

##@package convert_arcsec_pix
# 
#
# This package contains functions convert from arcsec units to pixels and vice-versa.


def convert_arcsec_pix(x,pix_scale):
	'''
	Convert the units from arcsec to pixels.
	
	Input:
	- x <float>: value in arcsec
	- pix_scale <float>: pixel scale (arcsec/pix)
	
	Output:
	- <float>: value in pixels
	'''
	return x / pix_scale


def convert_pix_arcsec(x,pix_scale):
	'''
	Convert the units from pixels to arcsec.
	
	Input:
	- x  <float>: value in pixel
	- pix_scale <float>: pixel scale (arcsec/pix)
	
	Output:
	- <float>: value in arcsec
	'''
	return x * pix_scale


def pix_2_arcsec(x_pix, y_pix, half_frame_size, dimpix):
    """
    Convert from pixel to arcsec coordinates. 

    The x and y coordinates in arcsec are calculated relative to the center of the image.
    
    Input:
     - x_pix          <ndarray> : 1d array (int) of the source x coordinates (in pixels)
     - y_pix          <ndarray> : 1d array (int) of the source y coordinates (in pixels)
     - half_frame_size  <float> : half the size of the FITS file (in arcsec)
     - dimpix           <float> : size of the pixel, in arcsec

    Output:
     - x_arcsec  <ndarray> : 1d array (int) of the source x coordinates (in arcsec)
     - y_arcsec  <ndarray> : 1d array (int) of the source y coordinates (in arcsec)

    """
    
    Npix = int( (2.*half_frame_size) / dimpix )
    x_arcsec = (2.*half_frame_size/(Npix-1))*x_pix - half_frame_size -2.*half_frame_size/(Npix - 1) # arcsec x coordinate
    y_arcsec = (2.*half_frame_size/(Npix-1))*y_pix - half_frame_size -2.*half_frame_size/(Npix - 1) # arcsec y coordinate

    return x_arcsec, y_arcsec 
