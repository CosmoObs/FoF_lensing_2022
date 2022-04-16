#!/usr/bin/env python
# ==================================
# Authors:
# Carlos Brandt - chbrandt@lncc.br
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Module to add noise to the image data """

##@package add_noise
#
# Module to add noise to the image data


import numpy as np
from numpy.random import poisson


def add_poisson_noise(img_array):

    """
    Add poissonian noise to given image array.

    Each nonzero pixel has its value updated from a Poisson distribution.    

    Input:
     - img_array <ndarray> : numpy.ndarray(ndim=2,dtype=float). Input data of an image
	
    Output:
     - <ndarray> : numpy.ndarray(ndim=2,dtype=float). The input array with Poisson noise added

    """

    ###########################################################################
    # Here I am treating a strange "behavior" (negative counts)
    #
    Arc_tmp = np.zeros( img_array.shape, dtype=float );
    non_null = np.where( img_array > 0.0 );
    Arc_tmp[ non_null ] = img_array[ non_null ];
    img_array = Arc_tmp.copy();
    del Arc_tmp;
    ###########################################################################


    # Re-sort pixel intensities of arc image in a poissonian fashion
    #
    img_array = poisson( img_array );


    return img_array;


