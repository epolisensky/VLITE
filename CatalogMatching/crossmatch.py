import numpy as np
from astropy.coordinates import SkyCoord


def deRuitermatch(src, catalog, bmaj, match_der=5.68, min_der=99999.9):
    """This matching function uses the unitless deRuiter
    radius to determine if there is a successful cross-match
    to a catalog source. A successful match must have a
    deRuiter radius < match_der, which defaults to 5.68, and
    have an angular distance < 0.5 * image beam width."""
    match = 0
    for catsrc in catalog:
        # If catalog is a list of catalog_source objects,
        # cast to a dictionary so it's attributes can be
        # accessed like the database row dictionary
        try:
            catsrc['dec']
        except TypeError:
            catsrc = catsrc.__dict__
        # Start by checking if within 15' in declination
        if quickcheck(src['dec'], catsrc['dec'], 0.25) == 1:
            # Calculate deRuiter radius
            der = deruiter(src['ra'], src['dec'], src['e_ra'], src['e_dec'],
                           catsrc['ra'], catsrc['dec'],
                           catsrc['e_ra'], catsrc['e_dec'])
            # Check if deRuiter radius < previous minimum
            if der < min_der:
                min_der = der
            # Check if deRuiter radius < required match limit
            if der < match_der:
                # Check if angular distance < 0.5*beam
                if angdist(src['ra'], src['dec'],
                           catsrc['ra'], catsrc['dec']) < (0.5 * bmaj):
                    match = 1
                    break
    return match, min_der, catsrc


def quickcheck(dec1, dec2, deglim):
    """Calculate the angular distance in declination 
    (in degrees) between two points (sources) and flag 
    if greater than the specified limit in degrees."""
    flag = 0
    if (abs(dec2 - dec1) < deglim):
      flag = 1
    return flag


def deruiter(ra1, dec1, e_ra1, e_dec1, ra2, dec2, e_ra2, e_dec2):
    """Calculate the unitless deRuiter radius between two 
    points (sources) using RA & Dec coordinates and their 
    errors in degrees."""
    r = np.sqrt(((((ra1 - ra2)**2.0) * (cosd((dec1 + dec2) / 2.0))**2.0) / \
                 (e_ra1**2.0 + e_ra2**2.0)) + (((dec1 - dec2)**2.0) / \
                                               (e_dec1**2.0 + e_dec2**2.0)))
    return r


def cosd(angle):
    """Returns the cosine of an angle given in degrees."""
    cosa = np.cos(np.radians(angle))
    return cosa


def angdist(ra1, dec1, ra2, dec2):
    """Return the angular distance (in degrees) 
    between two points (sources) given their RA & 
    Dec coordinates in degrees."""
    p1 = SkyCoord(ra1, dec1, unit="deg")
    p2 = SkyCoord(ra2, dec2, unit="deg")
    sep = p1.separation(p2)
    return sep.degree
