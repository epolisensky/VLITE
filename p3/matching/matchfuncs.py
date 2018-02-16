"""This script contains functions for positional cross-matching."""
import numpy as np
from astropy.coordinates import SkyCoord


def deRuitermatch(src, match_srcs, beam, match_der=6.44, min_der=99999.9):
    """This matching function uses the unitless de Ruiter
    radius to determine if there is a successful cross-match
    between two sources. A successful match must have a
    de Ruiter radius < match_der, which defaults to 6.44, and
    have an angular distance < 0.5 * image beam width in arcsec.
    
    Parameters
    ----------
    src : database.dbclasses.DetectedSource instance
        DetectedSource object to be matched.
    match_srcs : list
        List of either DetectedSource or CatalogSource
        objects which are candidates for matching to the
        given source.
    beam : float
        Semi-major or semi-minor axis of the image clean
        beam in arcseconds.
    match_der : float, optional
        The unitless minimum de Ruiter distance between two
        sources for them to be considered a match. Default
        is 6.44.
    min_der : float, optional
        Starting minimum de Ruiter distance. Default is 99999.9.

    Returns
    -------
    match : boolean
        Returns ``True`` if a match was successful for 
        the given source.
    match_src : object
        Returns the DetectedSource or CatalogSource object
        which was successfully matched to the given source.
        Otherwise, returns ``None``.
    min_der : float
        The minimum de Ruiter radius found for the given source.
    """
    match = False
    match_src = None
    for msrc in match_srcs:
       # Start by checking if within 15' in declination
        if quickcheck(src.dec, msrc.dec, 0.25):
            # Calculate de Ruiter radius
            der = deruiter(src.ra, src.dec, src.e_ra, src.e_dec,
                           msrc.ra, msrc.dec, msrc.e_ra, msrc.e_dec)
            # Check if de Ruiter radius < previous minimum
            if der < min_der:
                min_der = der
                match_src = msrc
            # Check if de Ruiter radius < required match limit
            if der < match_der:
                # Check if angular distance < 0.5*beam
                if angdist(src.ra, src.dec,
                           msrc.ra, msrc.dec) < (0.5 * beam):
                    match = True
                    break

    return match, match_src, min_der


def quickcheck(dec1, dec2, deglim):
    """Calculates the angular distance in degrees between
    the declinations of two points and sets the flag to
    ``True`` if less than the specified maximum separation
    limit in degrees.
    
    Parameters
    ----------
    dec1 : float
        Declination of the first source in degrees.
    dec2 : float
        Declination of the second source in degrees.
    deglim : float
        Maximum allowed separation between the two
        source's declinations, in degrees.

    Returns
    -------
    flag : boolean
        Returns ``True`` if the separation in declinations
        is less than the maximum allowed limit. Otherwise,
        returns ``False``.
    """
    flag = False
    if (abs(dec2 - dec1) < deglim):
      flag = True
      
    return flag


def deruiter(ra1, dec1, e_ra1, e_dec1, ra2, dec2, e_ra2, e_dec2):
    """Calculate the unitless de Ruiter radius between two 
    points using RA & Dec coordinates and their errors in degrees.
    
    Parameters
    ----------
    ra1 : float
        Right ascension in degrees of the primary source
        being matched.
    dec1 : float
        Declination in degrees of the primary source being
        matched.
    e_ra1 : float
        Error on the right ascension of the first source
        in degrees.
    e_dec1 : float
        Error on the declination of the first source in degrees.
    ra2 : float
        Right ascension in degrees of the second source attempting
        to be matched to the primary source.
    dec2 : float
        Declination in degrees of the second source attempting
        to be matched to the primary source.
    e_ra2 : float
        Error on the right ascension of the second source
        in degrees.
    e_dec2 : float
        Error on the declination of the second source in degrees.

    Returns
    -------
    r : float
        Calculated unitless de Ruiter distance between the two
        given sources.
    """
    r = np.sqrt(((((ra1 - ra2)**2.0) * (cosd((dec1 + dec2) / 2.0))**2.0) / \
                 (e_ra1**2.0 + e_ra2**2.0)) + (((dec1 - dec2)**2.0) / \
                                               (e_dec1**2.0 + e_dec2**2.0)))

    return r


def cosd(angle):
    """Returns the cosine of an angle given in degrees."""
    cosa = np.cos(np.radians(angle))
    
    return cosa


def angdist(ra1, dec1, ra2, dec2):
    """Returns the angular distance in arcsec between two points.

    Parameters
    ----------
    ra1 : float
        Right ascension of the first source in degrees.
    dec1 : float
        Declination of the first source in degrees.
    ra2 : float
        Right ascension of the second source in degrees.
    dec2 : float
        Declination of the second source in degrees.

    Returns
    -------
    sep : float
        The separation between the two points on the sky
        in arcseconds.
    """
    p1 = SkyCoord(ra1, dec1, unit="deg")
    p2 = SkyCoord(ra2, dec2, unit="deg")
    sep = p1.separation(p2).arcsecond
    
    return sep
