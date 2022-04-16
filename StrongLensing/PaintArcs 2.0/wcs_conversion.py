#!/usr/bin/python
# ====================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# ====================================================

"""Module to convert sky coordinates to pixel coordinates and vice-versa."""

##@package wcs_conversion
#
# This package contains functions to transform sky coordinates (RA and Dec) to pixel coordinates 
# (x and y) and vice-versa. It is a simple implementation of pywcs functions.
#
# 'radec2xy' transforms sky coordinates (RA and Dec) to pixel coordinates (x and y). 
#
# 'xy2radec' transforms pixel coordinates (x and y) to sky coordinates (RA and Dec).
#
# Both radec2xy and xy2radec have three mandatory arguments: the FITS image header and the coordinate 
# values (RA and Dec our x and y). RA and Dec values must be in degrees.
#
#

import os
import sys
import astropy.io.fits as pyfits
import numpy

try:
	import pywcs
except ImportError:
	import astropy.wcs as pywcs


def radec2xy(hdr,ra,dec):

	"""Transforms sky coordinates (RA and Dec) to pixel coordinates (x and y).
	
	Input:
	- hdr: FITS image header
	- ra <float> : Right ascension value in degrees
	- dec <float>: Declination value in degrees
	
	Output:
	- (x,y) <tuple>: pixel coordinates

	"""

	wcs = pywcs.WCS(hdr)
	
	skycrd = numpy.array([[ra,dec]])

	pixcrd = wcs.wcs_sky2pix(skycrd,1)

	x = pixcrd[0][0]

	y = pixcrd[0][1]

	return (x,y)
	
	


def xy2radec(hdr,x,y):

	""" Transforms pixel coordinates (x and y) to sky coordinates (RA and Dec).
	
	Input:
	- hdr: FITS image header
	- x <float>: x value
	- y <float>: y value
	
	Output:
	- (ra,dec) <tuple>: Right ascension and declination values in degrees 

	"""

	wcs = pywcs.WCS(hdr)

	pixcrd = numpy.array([[x,y]], numpy.float_)

	skycrd = wcs.wcs_pix2sky(pixcrd,1)

	ra = skycrd[0][0]

	dec = skycrd[0][1]

	return (ra,dec)	




