#!/usr/bin/env python

"""Module to deal with objects identification in segmented images"""

##@ segobjs
#
#
# This package contains functions to deal with image objects segmentation.
# Functions are designed primarily to work with Sextractor
#
# Executable package : No


import sys;

import astropy.io.fits as pyfits;
import numpy as np;


# =======================================================================
def centroid2id(segimg, centroid):
    """Determine ID of the objects which contains the given pixel.

    Input:
     - segimg : ndarray(ndim=2,dtype=int)
       ndarray with a segmented image (int), e.g, SE's SEGMENTATION 
     - centroid : (int,int)
       List of tuples with position to search for object (x0,y0)

    Output:
     - object ID : int
       List of IDs (int) for each given (x,y) points
       
       
    Example:
    >>> x = 1267
    >>> y = 851
    >>> centroid = zip(x,y)
    >>> objIDs = centroid2ID(SegImg_array,centroid)
    >>> objIDs
    2
    >>> 
    
    """


    # Is there an identified object on given (xo,yo) point? If yes, store its ID..
    #
    xo,yo = centroid;
    xo = int(float(xo));
    yo = int(float(yo));
    objid = segimg[yo,xo];

    return (objid);


#=============================================================================
def readout_ids(segimg):
    """ Read segmented image values as object IDs.
    
    Input:
     - segimg : ndarray
        Segmentation image array
    
    Output:
     -> list with object IDs : [int,]
        List with (typically integers) object IDs in 'segimg'
        Value "0" is taken off from output list ID
    
    """

    objIDs = list(set(seg_img.flatten()) - set([0]))

    return objIDs

#=============================================================================
def obj_id2index_array(segimg, objID):
    """ 
    Generate mask for each object in given image.
    
    ID in 'objIDs' is used to create a mask for each object (ID) in 'segimg'.
    
    Input:
     - segimg : ndarray(ndim=2,dtype=int)
        Image array with int numbers as object identifiers
     - objIDs : [int,]
        List with IDs for objects inside 'segimg'
        
    Output:
     -> index array (output from numpy.where())
        List of tuples with index arrays in it. The output is a list of "numpy.where()" arrays.
        
    """


    # For each object (id) scan the respective indices for image mask and 'cutout' params
    #
    id = float(objID);
    mask = np.where(segimg == int(id));
    
    return mask;

create_IDmask = obj_id2index_array;
