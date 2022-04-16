#!/usr/bin/env python
# ==================================
# Authors:
# Carlos Brandt - chbrandt@lncc.br
# ==================================

"""Package to make FITS image copy/cuts and update their header"""

##@package image
##@file imcp
#
#
# This package contains functions to cutout smaller regions
# from a original image, paste the image into a bigger array,
# as well as any other function for image copy manipulation.
#
# Executable package: YES
#
# To see the package help message, type:
#
# > python imcp.py --help
#
# 'imcp' has one mandatory argument: the (input) image FITS file.
# If none of the options are given, imcp outputs a file named
# 'output.fits' with 1/4th of the original image area, with the 
# same central point/pixel.
#
# > python imcp.py input_file.fits
#
# The options available make possible to define the central point
# where the cut will be centered, as well as the side lengths for
# output image. These values can be passed in 'pixel' or 'degrees'
# units.

import sys;
import logging;

import os
import astropy.io.fits as pyfits
import string
#import commands
import numpy as np;

import regexp;
import segmentation_ids;
import header_funcs as hf;
import wcs_conversion as wcsc;

class Interval:
    def __init__(self, begin = 0, end = 0):
        self.begin = begin
        self.end = end

def _cutshape(image, x_size, y_size, dimpix, unit):
    x_size = float(x_size)
    y_size = float(y_size)

    if x_size == 0:
        x_size = int(image.shape[0] / 2)
        logging.warning("'x_size' not given. Using half of image x side (%d)", x_size)

    if y_size == 0:
        y_size = int(image.shape[1] / 2)
        logging.warning("'y_size' not given. Using half of image y side (%d)", y_size)

    if unit == 'degrees':
        x_size = int(x_size / dimpix)
        y_size = int(y_size / dimpix)

    logging.debug("Output image shape: %s", (x_size, y_size))
    return x_size, y_size

def _cutcenter(image, hdr, xc, yc, coord_unit):
    if xc == 0 and yc == 0:
        xc = int(image.shape[0] / 2)
        yc = int(image.shape[1] / 2)
        logging.warning("No central coordinates were given for snapshot. Using image central pixel as the cut center.")
    else:
        if coord_unit == 'pixel':
            pass
        elif coord_unit == 'degrees':
            xc, yc = wcsc.radec2xy(hdr, xc, yc)
        else:
            logging.error("Central positioning is out of valid values (%d, %d)" % (xc, yc))

    return int(xc), int(yc)

def _image_patches(image, hdr, x, y, x_size, y_size):
    ox = Interval()
    ox.begin = int(max(0, x.begin))
    ox.end = int(min(image.shape[0], x.end))

    nx = Interval()
    nx.begin = int(abs(min(0, x.begin)))
    nx.end = int(x_size - (x.end - ox.end))

    oy = Interval()
    oy.begin = int(max(0, y.begin))
    oy.end = int(min(image.shape[0], y.end))

    ny = Interval()
    ny.begin = int(abs(min(0, y.begin)))
    ny.end = int(y_size - (y.end - oy.end))

    if hdr:
        hdr = hf.update_coordinates(hdr.copy(), x.begin, y.begin)
        hdr.update(NAXIS1 = x_size)
        hdr.update(NAXIS2 = y_size)

    return ox, oy, nx, ny

def cutout(image, hdr = None, coord_unit = 'pixel', xc = 0, yc = 0,
    size_unit = 'pixel', x_size = 0, y_size = 0, mask = None):
    """
    Produce a FITS image from a slice of another FITS image.

    cutout(ndarray, ...) -> (ndarray, header)

    An image mask can be given in a numpy.where like structure to denote a
    region of interest. The resulting image will have all pixels null except
    the ones listed in the 'mask' parameter.

    If 'xc = yc = 0', the center of the new image slice is the central pixel of
    the old one.

    If 'x_size = y_size = 0', the new image has half the width and height of
    the old one.


    Input:
     - image       numpy.ndarray : Image array (ndim = 2, dtype = float)
     - hdr         pyfits.header : Image FITS header
     - coord_unit            str : Position (xc, yc) units ('pixel', 'degrees')
     - xc                    int : Centroid (horizontal) position for the image slice
     - yc                    int : Centroid (vertical) position for the image slice
     - size_unit             str : Slice sizes (x_size, y_size) units ('pixel', 'degrees')
     - x_size                int : Horizontal slice size
     - y_size                int : Vertical slice size
     - mask   (ndarray, ndarray) : Tuple with index arrays (Output like numpy.where())

    Output:
     - (ndarray, header) : Image array of pixel intensities and (updated) headerx
     or
     - (False, False)    : On error
    
    ---
    """

    xc = float(xc)
    yc = float(yc)
    x_size = float(x_size)
    y_size = float(y_size)

    dimpix = hdr and hf.get_pixelscale(hdr) or 1
    x_size, y_size = _cutshape(image, x_size, y_size, dimpix, size_unit)

    logging.info("Pixel_scale: %s", dimpix)
    logging.info("Input image shape: %s", (image.shape))
    logging.info("Output image shape: %s", (x_size, y_size))

    if x_size == image.shape[0] and y_size == image.shape[1]:
        logging.warning("Requested output sizes are the same as input image. "
						"Returning image and header unchanged.")
        return image, hdr

    xc, yc = _cutcenter(image, hdr, xc, yc, coord_unit)
    if not xc:
        return False, False
    logging.info("Central point for output image: %s", (xc, yc))

    x = xc - int(x_size / 2)
    x = Interval(x, x + x_size)

    y = yc - int(y_size / 2)
    y = Interval(y, y + y_size)

    oldx, oldy, newx, newy = _image_patches(image, hdr, x, y, x_size, y_size)

    newimage = np.zeros((int(y_size), int(x_size)), dtype = image.dtype)
    newimage[newy.begin : newy.end, newx.begin : newx.end] = image[oldy.begin : oldy.end, oldx.begin : oldx.end]

    if hdr:
        hdr['XMIN'] = oldx.begin + 1
        hdr['YMIN'] = oldy.begin + 1

    #
    # If 'mask', maintain just "central" object on it..
    #
    if mask:
        ind_z = np.where(newimage == 0)
        mask = mask[0] - y.begin, mask[1] - x.begin

        zip_m = zip(mask[0], mask[1])
        zip_z = zip(ind_z[0], ind_z[1])

        L = list(set(zip_z) - set(zip_m))

        try:
            ind_0, ind_1 = zip(*L)
            indx = np.array(ind_0), np.array(ind_1)
            newimage[indx] = 0
        except:
            pass

    return newimage, hdr


def segstamp(segimg, id, objimg = None, hdr = None, increase = 0,
    relative_increase = False, connected = False, centered = True):
    """
    Identify objects on given images by their IDs and return object images

    segstamp( segimg, objID ... )

    By default, if 'objIDs' is not given, postamp will scan segmentation image 
    'seg_img' for the list of object ID numbers. If 'objIDs' is given, those IDs 
    will be used for object poststamps creation.

    'increase' and 'relative_increase' define whether the poststamps will have a 
    size different from object-only dimensions or just the object (pixels) itself, 
    and, also, if this increasement value is a multiplicative factor (True) or an 
    additive one (False).

    Since a list with object IDs is given, a list with arrays, with each IDentified
    object, is returned.

    Input:
     - segimg : numpy.ndarray(ndim=2,dtype=int)
        Segmented image (e.g, SEx's segmentation image)
        
     - objIDs : [int,]
        List with object IDs of interest in 'segimg'.

     - objimg : numpy.ndarray(ndim=2,dtype=float)
        Objects image (e.g, SEx's objects image)

     - hdr : FITS header instance
        FITS header to be updated for each poststamp
        
     - increase : float
        Value for poststamp resizing (> 0)
        
     - relative_increase : bool
        Is 'increase' a additive value (default,False) or multiplicative one(True)?
        

    Output:
     - (ndarray,header)  : Image (poststamp) array and corresponding header

    ---
    """

    mask = np.where(segimg == id)
    ypos = mask[0]
    xpos = mask[1]

    try:
        xmin = min(xpos)
        ymin = min(ypos)
        xmax = max(xpos)
        ymax = max(ypos)
    except ValueError:
        xmin = ymin = 0
        xmax = segimg.shape[1]
        ymax = segimg.shape[0]

    xsize = xmax - xmin + 1
    ysize = ymax - ymin + 1

    if centered:
        yc = ysize/2 + ymin
        xc = xsize/2 + xmin

        if increase != 0:
            if relative_increase:
                xsize = xsize * increase
                ysize = ysize * increase
            else:
                xsize = xsize + 2*increase
                ysize = ysize + 2*increase
    else:
        yc = ysize/2
        xc = xsize/2

        pixels = np.where(segimg >= 0)
        ysize = max(pixels[0])
        xsize = max(pixels[1])

    if connected:
      mask = separate_disconected(mask, high_area = True)

    if objimg.any():
        image = objimg
    else:
        image = segimg

    newimage, hdr = cutout(image, hdr, xc = int(xc), yc = int(yc),
	    x_size = int(xsize), y_size = int(ysize), mask = mask)

    return newimage, hdr

# ---
def select(segimg,objimg,objID):
    """
    Select pixels from a specific ID and outputs an array with the same input size
    
    Input:
     - segimg  ndarray : Segmented image array (ndim=2,dtype=int)
     - objimg  ndarray : Objects (observed) image array (ndim=2,dtype=float)
     - objID       int : Object ID to request output for
    
    Output:
     - outimg  ndarray : Same size (segimg) array with non-null pixels only for objID
    
    ---
    """

    outimg = np.zeros(segimg.shape)
    
    indxs = np.where(segimg==objID)
    if not len(indxs):
        return None
    
    outimg[indxs] = objimg[indxs]
    
    return outimg

# ---

def separate_disconected(ind, high_area=False):

    """
    From a list of points that belongs to objects it separates in groups of connected objects. If high_area=False return a list which first index represents each object , the second and third its coordinates. if high_area=True  returns the only a list with the coordinates of the object that has greater area.

   
    Input:
     - ind : array like
        the list of points where ind[0] is a list with the first coordinates (e.g. x) and ind[1] the second coordinates (e.g. y)
     - high_area : bool
        enables a criteria to return only the coordinates of the object that has greater area
   
    Output:
     - (nparray)  : if high_area=False  is list which first index represent the  object.  The second  and third index represents  the first and second coordinates. if high_area=True  is list which first and second index represents lists the first and second coordinates

    ---
    """

    p=[]
    for i in range(0,len(ind[0])):
        p_aux=[ind[0][i],ind[1][i]]
        p.append(p_aux)
   
  
    Objects=[]
    objects=[]
    while len(p)>0:
        p_test=[p[0]]
        while len(p_test)>0:
            p_center=p_test[0]
            p_neighbors=[[p_center[0],p_center[1]],[p_center[0]+1,p_center[1]],[p_center[0]+1,p_center[1]+1],[p_center[0]+1,p_center[1]-1],[p_center[0]-1,p_center[1]],[p_center[0]-1,p_center[1]+1],[p_center[0]-1,p_center[1]-1],[p_center[0],p_center[1]-1],[p_center[0],p_center[1]+1]]
            for i in range(0,len(p_neighbors)):
                if (p_neighbors[i] in p) and not(p_neighbors[i] in objects):
                    objects.append(p_neighbors[i])
                    p_test.append(p_neighbors[i])
                    p.remove(p_neighbors[i])  
            p_test.remove(p_test[0])
        Objects.append(objects)
        objects=[]

    if high_area==False:
        return Objects
    else:
        # criteria to select an interest object
   
        Area_max=len(Objects[0])
        id_max=0
        for i in range(0,len(Objects)):
            if len(Objects[i])>Area_max:
                Area_max=len(Objects[i])
                id_max=i
        ind_new=[[-100],[-1000]]
        for i in range(0,len(Objects[id_max])):
            ind_new[0].append(Objects[id_max][i][0])
            ind_new[1].append(Objects[id_max][i][1])
        ind_new[0].remove(-100)
        ind_new[1].remove(-1000)
        ind_n=[np.array(ind_new[0]),np.array(ind_new[1])]
        return ind_n          
