"""Functions for computing the primary beam correction."""

import numpy as np


def read_fitted_beam():
    """Reads the fitted beam file which contains
    the correction factor as a function of distance
    from the image center at every 10 arcseconds.
    
    Returns
    -------
    pbdir : dictionary
        Dictionary containing the primary beam
        correction factor ('power') as a function
        of the distance from the image center
        in degrees ('angle').
    """
    beamfile = '/home/erichards/work/data/FITBEAM_FINAL_DOMEGA.txt'
    with open(beamfile, 'r') as f:
        lines = f.readlines()

    pbdir = {'angle' : [], 'power' : [], 'domega' : []}
    for line in lines:
        data = line.split()
        pbdir['angle'].append(float(data[0]))
        pbdir['power'].append(float(data[1]))
        pbdir['domega'].append(float(data[2]))

    return pbdir


def find_nearest_pbcorr(angle):
    """Finds the primary beam correction factor
    corresponding to the angle value closest to 
    the given distance from the image center.

    Parameters
    ----------
    angle : float
        Distance from the image center in degrees.
    """
    pbdata = read_fitted_beam()
    idx = (np.abs(np.array(pbdata['angle']) - angle)).argmin()
    return pbdata['power'][idx]
