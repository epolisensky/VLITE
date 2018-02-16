"""Functions for computing the primary beam correction
and expected source counts.

"""
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
    #beamfile = '/home/vpipe/vlite-emily/data/FITBEAM_FINAL_DOMEGA.txt'
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

    Returns
    -------
    pbdata['power'][idx] : float
        Primary beam correction factor.
    """
    pbdata = read_fitted_beam()
    idx = (np.abs(np.array(pbdata['angle']) - angle)).argmin()
    
    return pbdata['power'][idx]


def expected_nsrc(rms, max_angle=1.5, sigma=5.):
    """Calculates expected number of sources that would
    have been detected within a maximum distance from the
    image center (max_angle) based on the detection threshold
    (sigma) in the image from logN-logS fit of WENSS catalog 
    scaled to 341 MHz from 325 MHz with alpha= -0.7, rms at 
    center of image [mJy], and fitted beam.

    Parameters
    ----------
    rms : float
        Estimate of the noise in the center of the image
        in mJy/beam.
    max_angle : float, optional
        Maximum angular distance in degrees from image center
        within which to count/estimate the number of sources.
        Default value is 1.5 degrees.
    sigma : float, optional
        The detection threshold used when source finding.
        Default value is 5.

    Returns
    -------
    nexp : int
        The expected number of sources within 'max_angle'
        given the detection threshold, 'sigma', and the
        image noise.
    """
    wenss_n0 = 10.637 # number per square degree
    wenss_S0 = 247.569 # mJy
    wenss_alpha1 = 1.846
    wenss_alpha2 = 0.332
    wenss_n = 0.497
    wenss_n2 = -1.0 / wenss_n
    # Read fitted beam file
    pbdata = read_fitted_beam()
    nexp = 0
    for i in range(len(pbdata['angle'])):
        if pbdata['angle'][i] > max_angle:
            break
        flux = ((sigma * rms) / pbdata['power'][i]) / wenss_S0
        nsrc = pbdata['domega'][i] * wenss_n0 * (
            (flux**(wenss_alpha1*wenss_n) \
             + flux**(wenss_alpha2*wenss_n))**wenss_n2)
        nexp += nsrc

    return nexp
