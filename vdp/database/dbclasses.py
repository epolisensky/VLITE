"""This module defines classes and their methods for
Image and DetectedSource objects.

"""
import sys
sys.path.insert(0,'/home/vpipe/VLITE/VLITE/vdp')
import os
import re
import logging
import numpy as np
import healpy as hp
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time
from astropy import wcs
import ephem
from database import pybdsfcat
from sourcefinding import beam_tools
from math import sqrt,ceil
from datetime import datetime

########################
# temporarily needed while USNO sites are down for modernization
#  expected completeion 30 Apr 2020
from astropy.utils import iers
#iers.conf.iers_auto_url = u'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
########################

RAD2DEG = 180.0/np.pi
DEG2RAD = np.pi/180.0

'''
VLALAT    = 34.078749   #VLA latitude  [deg]
VLALON    = -107.617728 #VLA longitude [deg]
VLAHEIGHT = 2124.0      #VLA height above sea level [m]
locVLA=EarthLocation(lat=VLALAT*u.deg, lon=VLALON*u.deg, height=VLAHEIGHT*u.m)
'''
locVLA = EarthLocation.of_site('vla')
VLALON = locVLA.lon.deg
VLALAT = locVLA.lat.deg

# create logger
dbclasses_logger = logging.getLogger('vdp.database.dbclasses')

# Define resolution dictionary
# Contains array configs, mjd ranges of the cycles, NRAO semesters and
#  beam solid angle normalization,
#  the acceptable range of image bmins for source association
res_dict = {
    'A': {'V10_1': {'mjd': [57182, 57212], 'bsanorm': 40.293, 'bmin': [0, 0], 'semester': ''},
          '1': {'mjd': [58179, 58280], 'bsanorm': 13.892750, 'bmin': [2.85, 4.27], 'semester': '2018A'},
          '2': {'mjd': [58697, 58778], 'bsanorm': 18.526056, 'bmin': [3.26, 4.90], 'semester': '2019A'},
          '3': {'mjd': [59190, 59285], 'bsanorm': 19.693765, 'bmin': [3.57, 5.36], 'semester': '2020B'},
          '4': {'mjd': [59648.1, 59766], 'bsanorm': 16.536189, 'bmin': [3.13, 4.70], 'semester': '2022A'},
          '5': {'mjd': [60129, 60220], 'bsanorm': 16.096211, 'bmin': [-1,-1], 'semester': '2023A'},
          '6': {'mjd': [60594, 99999], 'bsanorm': 17.839476, 'bmin': [-1,-1], 'semester': '2024B'}},
    'BnA': {'V10_1': {'mjd': [57154, 57182], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            'V10_2': {'mjd': [57639, 57668], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '1': {'mjd': [58148, 58179], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '2': {'mjd': [58660, 58697], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '3': {'mjd': [59142, 59190], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '4': {'mjd': [59610.4, 59648.1], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '5': {'mjd': [60094.46, 60129], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '6': {'mjd': [60570, 60594], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''}},
    'B': {'V10_1': {'mjd': [57050, 57154], 'bsanorm': 199.050, 'bmin': [0, 0], 'semester': ''},
          'V10_2': {'mjd': [57528, 57639], 'bsanorm': 313.233, 'bmin': [0, 0], 'semester': ''},
          '1': {'mjd': [57997, 58148], 'bsanorm': 203.477612, 'bmin': [10.77, 16.16], 'semester': '2017B'},
          '2': {'mjd': [58541, 58660], 'bsanorm': 155.607613, 'bmin': [10.18, 15.28], 'semester': '2019A'},
          '3': {'mjd': [59026, 59142], 'bsanorm': 231.048076, 'bmin': [11.27, 16.91], 'semester': '2020A'},
          '4': {'mjd': [59475.4, 59610.4], 'bsanorm': 211.639557, 'bmin': [10.33, 15.50], 'semester': '2021B'},
          '5.1': {'mjd': [59957, 59993], 'bsanorm': 175.283348, 'bmin': [-1,-1], 'semester': '2023A'},
          '5.2': {'mjd': [59993, 60094.46], 'bsanorm': 209.598754, 'bmin': [-1,-1], 'semester': '2023A'},
          '6': {'mjd': [60431, 60570], 'bsanorm': 212.462536, 'bmin': [-1,-1], 'semester': '2024A'}},
    'CnB': {'V10_1': {'mjd': [57029, 57050], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            'V10_2': {'mjd': [57526, 57528], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '1': {'mjd': [57994, 57997], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '2': {'mjd': [58519, 58541], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '3': {'mjd': [59008, 59026], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '4': {'mjd': [59471, 59475.4], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '5': {'mjd': [59948, 59957], 'bsanorm': -1, 'bmin': [-1,-1], 'semester': ''},
            '6': {'mjd': [60423, 60431], 'bsanorm': -1, 'bmin': [-1,-1], 'semester': ''}},
    'C': {'V10_1': {'mjd': [56986, 57029], 'bsanorm': 1910.885, 'bmin': [32.2, 65.91], 'semester': ''},
          'V10_2': {'mjd': [57423.6, 57526], 'bsanorm': 1967.009, 'bmin': [32.2, 65.91], 'semester': ''},
          'V10_3': {'mjd': [57898, 57954], 'bsanorm': 2134.905, 'bmin': [32.2, 65.91], 'semester': ''},
          '1': {'mjd': [57954, 57994], 'bsanorm': 2686.817993, 'bmin': [43.94, 65.91], 'semester': '2017A'},
          '2': {'mjd': [58441, 58519], 'bsanorm': 1744.137879, 'bmin': [32.20, 48.30], 'semester': '2018B'},
          '3': {'mjd': [58885, 59008], 'bsanorm': 2393.579216, 'bmin': [39.69, 59.54], 'semester': '2020A'},
          '4': {'mjd': [59419, 59471], 'bsanorm': 1730.897165, 'bmin': [33.29, 49.93], 'semester': '2021A'},
          '5': {'mjd': [59853, 59948], 'bsanorm': 2131.092858, 'bmin': [-1,-1], 'semester': '2022B'},
          '6': {'mjd': [60334, 60423], 'bsanorm': 1877.441356, 'bmin': [-1,-1], 'semester': '2024A'}},
    'DnC': {'V10_2': {'mjd': [57393, 57423.6], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            'V10_3': {'mjd': [57890, 57898], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '2': {'mjd': [58437, 58441], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '3': {'mjd': [58876, 58885], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '4': {'mjd': [59367, 59374], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '5': {'mjd': [59848, 59853], 'bsanorm': -1, 'bmin': [-1,-1], 'semester': ''},
            '6': {'mjd': [60324, 60334], 'bsanorm': -1, 'bmin': [-1,-1], 'semester': ''}},
    'D': {'V10_2': {'mjd': [57308, 57393], 'bsanorm': 13061.635, 'bmin': [0, 0], 'semester': ''},
          'V10_3': {'mjd': [57794, 57890], 'bsanorm': 16572.679, 'bmin': [0, 0], 'semester': ''},
          '2': {'mjd': [58360, 58437], 'bsanorm': 20271.136987, 'bmin': [117.58, 176.37], 'semester': '2018A'},
          '3': {'mjd': [58802, 58876], 'bsanorm': 24069.990270, 'bmin': [133.07, 199.60], 'semester': '2019B'},
          '4': {'mjd': [59295, 59367], 'bsanorm': 22982.854365, 'bmin': [127.37, 191.05], 'semester': '2021A'},
          '5': {'mjd': [59785, 59848], 'bsanorm': 21899.169321, 'bmin': [-1,-1], 'semester': '2022A'},
          '6': {'mjd': [60237, 60324], 'bsanorm': -1, 'bmin': [-1,-1], 'semester': '2023B'}}, #VLITE OFF!
    'B+': {'V10_1': {'mjd': [57212, 57295], 'bsanorm': 99.404, 'bmin': [0, 0], 'semester': ''},
           'V10_2': {'mjd': [57668, 57778], 'bsanorm': 60.267, 'bmin': [0, 0], 'semester': ''}},
    'A-D': {'V10_1-2': {'mjd': [57295, 57308], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            'V10_2-3': {'mjd': [57778, 57794], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '1-2': {'mjd': [58280, 58360], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '2-3': {'mjd': [58778, 58802], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '3-4': {'mjd': [59285, 59295], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '4-5': {'mjd': [59766, 59785], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''},
            '5-6': {'mjd': [60220, 60237], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''}},
    'Unk': {'0': {'mjd': [99998, 99999], 'bsanorm': -1, 'bmin': [0, 0], 'semester': ''}}
}

compactness_dict = {
    'A' : {'1' : {'a0': 1.031571, 'a1': 4.266747e-03, 'c0': 9.784101, 'c1': -1.666293},
           '2' : {'a0': 1.018140, 'a1': 2.699981e-03, 'c0': 6.541704, 'c1': -1.538327},
           '3' : {'a0': 1.019842, 'a1': 2.439492e-03, 'c0': 8.346418, 'c1': -1.556666},
           '4' : {'a0': 1.037543, 'a1': 5.210229e-03, 'c0': 14.28838, 'c1': -1.801248},
           '5' : {'a0': 1.057399, 'a1': 9.871675e-03, 'c0': 10.88056, 'c1': -1.763789},
           '6' : {'a0': 0, 'a1': 0, 'c0': 0, 'c1': -0},
           'V10_1' : {'a0': 1.0, 'a1': 2.831733e-04, 'c0': 3.346050, 'c1': -1.001673}},
    'B' : {'1' : {'a0': 1.021638, 'a1': 1.284321e-03, 'c0': 6.533840, 'c1': -1.617770},
           '2' : {'a0': 1.025000, 'a1': 1.960815e-03, 'c0': 5.327601, 'c1': -1.557729},
           '3' : {'a0': 1.013762, 'a1': 1.104968e-03, 'c0': 4.491889, 'c1': -1.451631},
           '4' : {'a0': 1.017186, 'a1': 1.274251e-03, 'c0': 4.754481, 'c1': -1.506705},
           '5.1' : {'a0': 1.005072, 'a1': 7.54856e-04, 'c0': 3.301659, 'c1': -1.341046},
           '5.2' : {'a0': 1.010506, 'a1': 9.91735e-04, 'c0': 3.481263, 'c1': -1.393633},
           '6' : {'a0': 1.002334, 'a1': 9.046608e-04, 'c0': 2.150379, 'c1': -1.236555},
           'V10_1' : {'a0': 1.007135, 'a1': 2.925095e-03, 'c0': 4.855276, 'c1': -1.499755},
           'V10_2' : {'a0': 1.030928, 'a1': 4.127990e-03, 'c0': 6.398501, 'c1': -1.633580}},
    'C' : {'1' : {'a0': 1.005082, 'a1': 8.677483e-04, 'c0': 3.476722, 'c1': -1.549343},
           '2' : {'a0': 1.010657, 'a1': 1.655885e-03, 'c0': 2.591503, 'c1': -1.365536},
           '3' : {'a0': 1.013413, 'a1': 1.218773e-03, 'c0': 2.570887, 'c1': -1.465882},
           '4' : {'a0': 1.021204, 'a1': 1.598263e-03, 'c0': 2.603568, 'c1': -1.499281},
           '5' : {'a0': 1.020900, 'a1': 1.759167e-03, 'c0': 2.626716, 'c1': -1.544575},
           '6' : {'a0': 1.033150, 'a1': 2.817737e-03, 'c0': 2.897327, 'c1': -1.553051},
           'V10_1' : {'a0': 1.006811, 'a1': 1.476478e-03, 'c0': 2.679202, 'c1': -1.463098},
           'V10_2' : {'a0': 1.012908, 'a1': 1.788574e-03, 'c0': 1.606308, 'c1': -1.148350},
           'V10_3' : {'a0': 1.003027, 'a1': 1.152652e-03, 'c0': 2.942781, 'c1': -1.387235}},
    'D' : {'2' : {'a0': 1.029931, 'a1': 3.368299e-03, 'c0': 3.831083, 'c1': -1.696777},
           '3' : {'a0': 1.011935, 'a1': 1.780312e-03, 'c0': 2.105211, 'c1': -1.401127},
           '4' : {'a0': 1.016972, 'a1': 1.515647e-03, 'c0': 2.492477, 'c1': -1.549401},
           '5' : {'a0': 1.020662, 'a1': 2.635404e-03, 'c0': 3.467282, 'c1': -1.627111},
           '6' : {'a0': 0, 'a1': 0, 'c0': 0, 'c1': -0},
           'V10_2' : {'a0': 0.0, 'a1': 0.0, 'c0': 0.0, 'c1': -0.0},
           'V10_3' : {'a0': 0.0, 'a1': 0.0, 'c0': 0.0, 'c1': -0.0}},
    'B+' : {'V10_1' : {'a0': 1.027355, 'a1': 2.882720e-03, 'c0': 12.065286, 'c1': -1.783792},
           'V10_2' : {'a0': 1.055854, 'a1': 8.093716e-03, 'c0': 15.708964, 'c1': -1.867890}},
    'VCSS' : {'1.1' : {'a0': 1.010937, 'a1': 7.324113e-04, 'c0': 2.917032, 'c1': -1.496858},
              '1.2' : {'a0': 1.008186, 'a1': 6.611323e-04, 'c0': 3.978883, 'c1': -1.532021},
              '2.1' : {'a0': 1.020981, 'a1': 8.210672e-04, 'c0': 7.564879, 'c1': -1.807970},
              '2.2' : {'a0': 1.018776, 'a1': 9.647108e-04, 'c0': 2.845967, 'c1': -1.668669},
              '3.1' : {'a0': 1.016622, 'a1': 7.122767e-04, 'c0': 1.295902, 'c1': -1.332420},
              '3.2' : {'a0': 0.0, 'a1': 0.0, 'c0': 0.0, 'c1': -0.0}},
    }

class ImgMjd(object):
    """A class to hold the name & mjd of a FITS image.
       Used for sorting images by mjd before processing fully.
    Parameters
    ----------
    image : str
        Directory path to the FITS image file location.

    Attributes
    ----------
    filename : str
        Full directory path for the image file.
    mjdtime : float
        Modified Julian Date (days since 0h Nov 17, 1858).
    """
    def __init__(self):
        self.filename = None
        self.mjdtime = None

    def init_image(self, image):
        self.filename = image
        try:
            hdu = fits.open(self.filename, mode='readonly')
            hdr = hdu[0].header
            obs_date = hdr['DATE-OBS']
            if len(obs_date) > 10:
                obs_date = hdr['DATE-OBS'][:10]
            try:
                if obs_date is not None:
                    date = obs_date.split('-')
                    try:
                        self.mjdtime = Time(datetime(int(date[0]), int(date[1]), int(
                            date[2]))).mjd + float(hdr['STARTIME'])  # day
                    except:
                        self.mjdtime = Time(datetime(int(date[0]), int(date[1]), int(
                            date[2]))).mjd  # day
            except KeyError:
                self.mjdtime = None
            hdu.close()
        except:
            dbclasses_logger.error('\nERROR: Problem reading image.')

    
            

class Image(object):
    """A class to hold information about the FITS image file.

    Parameters
    ----------
    image : str
        Directory path to the FITS image file location.

    Attributes
    ----------
    id : int
        Numerical index. Incremented by PostgreSQL upon insertion
        into the database **image** table.
    filename : str
        Full directory path for the image file.
    wcsobj : object
        image header converted to WCS object
    imsize : str
        Image size in pixels -- (``NAXIS1``, ``NAXIS2``).
    obs_ra : float
        Right ascension of image center (degrees). Keyword: CRVAL1
    obs_dec : float
        Declination of image center (degrees). Keyword: CRVAL2
    point_ra : float
        Right ascension of VLA pointing center (degrees). Keyword: OBSRA
    point_dec : float
        Declination of VLA pointing center (degrees). Keyword: OBSDEC
    glon : float
        Galactic longitude of image pointing center (degrees).
    glat : float
        Galactic latitude of image pointing center (degrees).
    az_star : float
        Azimuth of image pointing center at start time from header (degrees).
    el_star : float
        Elevation of image pointing center at start time from header (degrees).
    pa_star : float
        Parallactic angle of image pointing center at start time from header (degrees).
    az_end : float
        Azimuth of image pointing center at end time from header (degrees).
    el_end : float
        Elevation of image pointing center at end time from header (degrees).
    pa_end : float
        Parallactic angle of image pointing center at end time from header (degrees).
    az_i : float
        Azimuth of image pointing center at mjdtime from astropy (degrees).
    alt_i : float
        Altitude of image pointing center at mjdtime from astropy (degrees).
    parang_i : float
        Parallactic angle of image pointing center at mjdtime from astropy (degrees).
    hrang_i : float
        Hour angle at mjdtime (hrs).
    az_f : float
        Azimuth of image pointing center at mjdtime+duration from astropy (degrees).
    alt_f : float
        Altitude of image pointing center at mjdtime+duration from astropy (degrees).
    parang_f : float
        Parallactic angle of image pointing center at mjdtime+duration from astropy (degrees).
    hrang_f : float
        Hour angle at mjdtime+duration (hrs).
    lst_i : float
        Local Sidereal Time of image at mjdtime (hours) 
    lst_f : float
        Local Sidereal Time at mjdtime+duration (hours)
    pixel_scale : float
        Spatial conversion from pixel to angular size on the sky
        (arcseconds / pixel.)
    obj : str
        Name of observed object in the header.
    obs_date : date
        Date when observations were taken.
    map_date : date
        Date when data was imaged.
    obs_freq : float
        Rest frequency of the observations (MHz).
    pri_freq : float
        Frequency of the simultaneously acquired primary
        band data (GHz).
    priband : str
        String version of pri_freq.
    pbkey : str
        Primary beam key for image (depends on priband).
    bmaj : float
        Size of the beam major axis FWHM (arcsec).
    bmin : float
        Size of the beam minor axis FWHM (arcsec).
    bpa : float
        Beam position angle (degrees).
    noise : float
        Noise estimate from the image center (mJy/beam).
    peak : float
        Peak brightness in the image (mJy/beam).
    config : str
        Very Large Array configuration.
    cycle : str
        VLITE cycle corresponding to config
    semester: str
        NRAO semester corresponding to config
    nvis : int
        Number of visibilities in the data after calibration.
    niter : int
        Number of CLEAN iterations
    mjdtime : float
        Modified Julian Date (days since 0h Nov 17, 1858).
    tau_time : float
        Integration time on source (seconds).
    duration : float
        Total time spanned by the observations (seconds).
    nsrc : int
        Number of sources extracted during source finding.
    nclean : int
        Number of sources that were CLEANed
    rms_box : str
        PyBDSF parameter used during source finding.
        From PyBDSF docs: "The first integer, boxsize, is the 
        size of the 2-D sliding box [in pixels] for calculating 
        the rms and mean over the entire image. The second, stepsize, 
        is the number of pixels by which this box is moved for the next 
        measurement."
    error_id : int
        Numerical code assigned when an image fails a quality
        check. Each number has a corresponding explanation in
        the database **error** table.
    stage : int
        Highest processing stage completed. 1 = read image,
        2 = source finding, 3 = source association, 4 = catalog matching.
    radius : float
        Radius defining the circular field-of-view in which
        sources were extracted (degrees).
    nearest_problem : str
        Name of the closest (in angular separation) source which is known
        to cause difficulties in imaging.
    separation : float
        Angular separation between the nearest problem source
        and the image pointing center (degrees).
    pri_cals : json
        Primary calibrators used. Read from history extension
    cc : list of (float, float)
        (RA, Dec) list of clean components. Read from CC extension
    ass_flag : boolean
        If True image meets requirments for source association
    nsn : int
        Number of SN tables in UVOUT. 3 = amp & phase self cal applied to image 
    tsky : float
        Brightness temp of GSM map at image center
    square : str
        VCSS Square. 
    xref : float
        X reference pixel. Keyword CRPIX1
    yref : float
        y reference pixel. Keyword CRPIX2
    xpoint : float
        X reference pixel of VLA pointing center
    ypoint : float
        y reference pixel of VLA pointing center
    naxis1 : int
        Number of x-dimension pixels. Keyword NAXIS1
    naxis2 : int
        Number of y-dimension pixels. Keyword NAXIS2
    bmimg : numpy ndarray
        primary beam image for the image
    vcss : boolean
        True if image is a VCSS snapshot (used for beam image calculation)
    mosaic : boolean
        True if image is a VCSS mosaic
    sunsep : float
        Angular separation between the Sun and the image pointing 
        center at mjdtime (degrees).
    pb_flag : boolean
        True if able to calculate primary beam image
    ninterval : integer
        Number of time intervals in NX table
    max_dt : float
        Max time interval in NX table
    nvisnx : integer
        Number of visibilities from NX table. Before self-cal
    nbeam : integer
        Number of beams in beam image calculation 
    pbtimes : list of nbeam floats
        Times for beam image calculations (days)
    pbparangs : list of nbeam floats
        Parallactic angles of beam image calculation (deg)
    pbweights : list of nbeam floats
        Weights to apply to beams calculated at pbtimes
    pbza : list of nbeam floats
        Zenith angles of beam image calculation (deg)
    smeartime : float
        From setup, smear time of pri beam image calculation (s)
    startime : float
        Header keyword STARTIME [days]
    """
    # A class variable to count the number of images
    num_images = 0

    def __init__(self):
        self.id = None
        self.filename = None
        self.wcsobj = None
        self.imsize = None
        self.obs_ra = None
        self.obs_dec = None
        self.point_ra = None
        self.point_dec = None
        self.glon = None
        self.glat = None
        self.az_star = None
        self.el_star = None
        self.pa_star = None
        self.az_end = None
        self.el_end = None
        self.pa_end = None
        self.az_i = None
        self.alt_i = None
        self.parang_i = None
        self.az_f = None
        self.alt_f = None
        self.parang_f = None
        self.lst_i = None
        self.lst_f = None
        self.pixel_scale = None
        self.obj = None
        self.obs_date = None
        self.map_date = None
        self.obs_freq = None
        self.pri_freq = None
        self.bmaj = None
        self.bmin = None
        self.bpa = None
        self.noise = None
        self.peak = None
        self.config = None
        self.cycle = None
        self.semester = None
        self.nvis = None
        self.niter = None
        self.mjdtime = None
        self.tau_time = None
        self.duration = None
        self.nsrc = None
        self.nclean = None
        self.rms_box = None
        self.error_id = None
        self.stage = 1
        self.radius = None
        self.nearest_problem = None
        self.separation = None
        self.pri_cals = None
        self.cc = None
        self.ass_flag = None
        self.nsn = None
        self.tsky = None
        self.square = None
        self.xref = None
        self.yref = None
        self.xpoint = None
        self.ypoint = None
        self.naxis1 = None
        self.naxis2 = None
        self.bmimg = None
        self.vcss = None
        self.mosaic = None
        self.sunsep = None
        self.hrang_i = None
        self.hrang_f = None
        self.pbkey = None
        self.pb_flag = None
        self.ninterval = None
        self.max_dt = None
        self.nvisnx = None
        self.nbeam = None
        self.pbtimes = None
        self.pbweights = None
        self.pbparangs = None
        self.pbza = None
        self.smeartime = None
        self.startime = None

    def process_image(self, image):
        self.filename = image

    @classmethod
    def image_count(cls):
        """Returns the number of Image objects initialized."""
        return cls.num_images

    def read(self):
        """Reads FITS image data and header."""
        try:
            hdu = fits.open(self.filename, mode='readonly')
            hdr = hdu[0].header
            return hdu, hdr
        except:
            dbclasses_logger.error('\nERROR: Problem reading image.')

    # necessary if keywords needed by PyBDSF are missing or
    #  incorrectly labeled (CTYPE3)
    def write(self, hdu, header):
        """Updates FITS image header"""
        dbclasses_logger.info('Updating header of {}'.format(self.filename))
        #open in update mode
        hdu2 = fits.open(self.filename, mode='update')
        hdr2 = hdu2[0].header
        hdr2['CTYPE3'] = header['CTYPE3']
        hdr2['BMAJ'] = header['BMAJ']
        hdr2['BMIN'] = header['BMIN']
        hdr2['BPA'] = header['BPA']
        hdu2.flush()
        hdu2.close()

    def header_attrs(self, hdr):
        """Extracts all keywords of interest from the FITS image
        header and stores their values as attributes of the
        initialized Image object. If a header keyword is missing
        from the image metadata, then that attribute value is
        set to ``None``. The error_id attribute is also set to 1
        if the missing header keyword is deemed important enough.        

        """
        # first check if VCSS snapshot
        self.vcss = False
        if 'VCSS' in self.filename or 'vcss' in self.filename:
                if self.filename.endswith('IPln1.fits'):
                    self.vcss = True # image is a VCSS snapshot not mosaic
                    self.pb_flag = True

        # check if VCSS mosaic
        self.mosaic = False
        if self.filename.endswith('IMSC.fits'):
            self.mosaic = True # image is a VCSS mosaic
            self.pb_flag = False
        
        # start with priband
        try:
            hdrpriband = hdr['PRIBAND']
            if float(hdrpriband[:-3]) > 100:
                self.pri_freq = float(hdrpriband[:-3])/1000  # GHz
                self.priband = '0.3'
                self.pbkey = 'P'
            else:
                self.pri_freq = float(hdrpriband[:-3])  # GHz
                if self.pri_freq > 40:
                    self.priband = '45'
                    self.pbkey = 'KQ'
                elif self.pri_freq > 30:
                    self.priband = '33'
                    self.pbkey = 'KQ'
                elif self.pri_freq > 20:
                    self.priband = '22'
                    self.pbkey = 'KQ'
                elif self.pri_freq > 14:
                    self.priband = '15'
                    self.pbkey = 'KQ'
                elif self.pri_freq > 9:
                    self.priband = '10'
                    self.pbkey = 'CX'
                elif self.pri_freq > 5:
                    self.priband = '6'
                    self.pbkey = 'CX'
                elif self.pri_freq > 2:
                    self.priband = '3'
                    self.pbkey = 'S'
                elif self.pri_freq > 1:
                    self.priband = '1.5'
                    self.pbkey = 'L'
                else:
                    self.priband = None
                    self.pbkey = None
        except:
            # VCSS mosaics & snapshots
            if self.mosaic or self.vcss:
                self.pri_freq = 3
                self.priband = '3'
                self.pbkey = 'S'
            else:
                pri = re.findall('\/([0-9.]+[A-Z])', self.filename)
                try:
                    if pri[0][-1:] == 'M':
                        self.pri_freq = float(pri[0][:-1]) / 1000.  # GHz
                        self.priband = '0.3'
                        self.pbkey = 'P'
                    else:
                        self.pri_freq = float(pri[0][:-1])  # GHz
                        if self.pri_freq > 40:
                            self.priband = '45'
                            self.pbkey = 'KQ'
                        elif self.pri_freq > 30:
                            self.priband = '33'
                            self.pbkey = 'KQ'
                        elif self.pri_freq > 20:
                            self.priband = '22'
                            self.pbkey = 'KQ'
                        elif self.pri_freq > 14:
                            self.priband = '15'
                            self.pbkey = 'KQ'
                        elif self.pri_freq > 9:
                            self.priband = '10'
                            self.pbkey = 'CX'
                        elif self.pri_freq > 5:
                            self.priband = '6'
                            self.pbkey = 'CX'
                        elif self.pri_freq > 2:
                            self.priband = '3'
                            self.pbkey = 'S'
                        elif self.pri_freq > 1:
                            self.priband = '1.5'
                            self.pbkey = 'L'
                        else:
                            self.priband = None
                            self.pbkey = None
                    # project code 13B-266 is at 1.5 GHz
                    if self.pri_freq == 13:
                        self.pri_freq = 1.5
                        self.priband = '1.5'
                        self.pbkey = 'L'
                except IndexError:
                    self.pri_freq = None
                    self.priband = None
                    self.pbkey = None
        ###
        try:
            self.naxis1 = hdr['NAXIS1']
            self.naxis2 = hdr['NAXIS2']
            self.imsize = str((self.naxis1, self.naxis2))  # pixels
        except KeyError:
            self.naxis1 = None
            self.naxis2 = None
            self.imsize = None
            self.error_id = 1

        self.wcsobj = wcs.WCS(hdr).celestial #keeps only RA, Dec axes

        try:
            self.obs_ra = hdr['CRVAL1']  # deg
        except:
            self.obs_ra = None
            self.error_id = 1
        try:
            self.point_ra = hdr['OBSRA']  # deg
        except KeyError:
            self.point_ra = None
            self.error_id = 1
        try:
            self.obs_dec = hdr['CRVAL2']  # deg
        except KeyError:
            self.obs_dec = None
            self.error_id = 1
        try:
            self.point_dec = hdr['OBSDEC']  # deg
        except KeyError:
            self.point_dec = None
            self.error_id = 1

        self.xpoint, self.ypoint = beam_tools.EQtoPIX(self.point_ra, self.point_dec, self.wcsobj)
        
        try:
            self.pixel_scale = abs(hdr['CDELT1']) * 3600.  # arcsec/pixel
        except KeyError:
            try:
                self.pixel_scale = abs(hdr['CDELT2']) * 3600.  # arcsec/pixel
            except KeyError:
                self.pixel_scale = None
                self.error_id = 1
        try:
            self.duration = hdr['DURATION']  # sec
        except KeyError:
            self.duration = None
            self.error_id = 1
        try:
            self.xref = hdr['CRPIX1']
        except KeyError:
            self.xref = None
            self.error_id = 1
        try:
            self.yref = hdr['CRPIX2']
        except KeyError:
            self.yref = None
            self.error_id = 1    
        try:
            self.obj = hdr['OBJECT']
        except KeyError:
            self.obj = None
        try:
            self.obs_date = hdr['DATE-OBS']
            if len(self.obs_date) > 10:
                self.obs_date = hdr['DATE-OBS'][:10]
        except KeyError:
            self.obs_date = None
            self.error_id = 1
        try:
            self.map_date = hdr['DATE-MAP']
        except KeyError:
            self.map_date = None
        try:
            self.obs_freq = hdr['RESTFREQ'] / 10**6.  # MHz
        except KeyError:
            try:
                if hdr['CTYPE3'] == 'FREQ' or hdr['CTYPE3'] == 'SPECLNMF':
                    self.obs_freq = hdr['CRVAL3'] / 10**6.  # MHz
                else:
                    self.obs_freq = hdr['CRVAL4'] / 10**6.  # MHz
            except KeyError:
                self.obs_freq = None
                self.error_id = 1
        try:
            self.bmaj = hdr['BMAJ'] * 3600.  # arcsec
            self.bmin = hdr['BMIN'] * 3600.  # arcsec
            self.bpa = hdr['BPA']  # deg
        except KeyError:
            try:
                self.bmaj = hdr['CLEANBMJ'] * 3600.  # arcsec
                self.bmin = hdr['CLEANBMN'] * 3600.  # arcsec
                self.bpa = hdr['CLEANBPA']  # deg
            except KeyError:
                try:
                    # Search for beam params in AIPS history
                    hl = list(hdr['HISTORY'])
                    for line in hl:
                        x = re.findall('BMAJ=\s+([0-9]\S+)', line)
                        y = re.findall('BMIN=\s+([0-9]\S+)', line)
                        z = re.findall('BPA=\s+([0-9]\S+)', line)
                        if len(x) > 0:
                            self.bmaj = float(x[0]) * 3600.  # arcsec
                        if len(y) > 0:
                            self.bmin = float(y[0]) * 3600.  # arcsec
                        if len(z) > 0:
                            self.bpa = float(z[0])  # deg
                except KeyError:
                    self.bmaj = None
                    self.bmin = None
                    self.bpa = None
                    self.error_id = 1
        try:
            self.noise = hdr['ACTNOISE'] * 1000.  # mJy/beam
        except KeyError:
            self.noise = None
            self.error_id = 1
        try:
            self.peak = hdr['PEAK'] * 1000.  # mJy/beam
        except KeyError:
            try:
                self.peak = hdr['DATAMAX'] * 1000.  # mJy/beam
            except KeyError:
                self.peak = None
        try:
            self.nvis = hdr['NVIS']
        except KeyError:
            self.nvis = None
            self.error_id = 1
        try:
            self.niter = hdr['CLEANNIT']
        except KeyError:
            try:
                self.niter = hdr['NITER']
            except KeyError:
                try:
                    hl = list(hdr['HISTORY'])
                    for line in hl:
                        x = re.findall('NITER=\s+([0-9]\S+)', line)
                        if len(x) > 0:
                            self.niter = int(x[0])
                except KeyError:
                    self.niter = None
                    self.error_id = 1
        try:
            if self.obs_date is not None:
                date = self.obs_date.split('-')
                try:
                    self.mjdtime = Time(datetime(int(date[0]), int(date[1]), int(
                        date[2]))).mjd + float(hdr['STARTIME'])  # day
                    self.startime = float(hdr['STARTIME'])
                except:
                    self.mjdtime = Time(datetime(int(date[0]), int(date[1]), int(
                        date[2]))).mjd  # day
                t = Time(self.mjdtime, format='mjd')
                self.lst_i = t.sidereal_time('apparent', VLALON).hour  # hrs
                if self.point_ra is not None and self.point_dec is not None:
                    coord = SkyCoord(self.point_ra, self.point_dec,
                                     unit='deg', frame='fk5')
                    altaz = coord.transform_to(
                        AltAz(obstime=t, location=locVLA))
                    self.az_i = altaz.az.deg
                    self.alt_i = altaz.alt.deg
                    hrang = (15*self.lst_i) - self.point_ra  # deg
                    if hrang > 180:
                        hrang -= 360
                    if hrang < -180:
                        hrang += 360
                    self.hrang_i = hrang/15 # hrs
                    tmp1 = np.sin(hrang*DEG2RAD)
                    tmp2 = np.tan(VLALAT*DEG2RAD)*np.cos(self.point_dec*DEG2RAD) - \
                        np.sin(self.point_dec*DEG2RAD)*np.cos(hrang*DEG2RAD)
                    self.parang_i = np.arctan2(tmp1, tmp2)*RAD2DEG
                    if self.duration is not None:
                        t_end = Time(self.mjdtime+(self.duration /
                                                   86400.), format='mjd')  # end time
                        self.lst_f = t_end.sidereal_time(
                            'apparent', VLALON).hour  # hrs
                        altaz = coord.transform_to(
                            AltAz(obstime=t_end, location=locVLA))
                        self.az_f = altaz.az.deg
                        self.alt_f = altaz.alt.deg
                        hrang = (15*self.lst_f) - self.point_ra  # deg
                        if hrang > 180:
                            hrang -= 360
                        if hrang < -180:
                            hrang += 360
                        self.hrang_f = hrang/15 # hrs
                        tmp1 = np.sin(hrang*DEG2RAD)
                        tmp2 = np.tan(VLALAT*DEG2RAD)*np.cos(self.point_dec*DEG2RAD) - \
                            np.sin(self.point_dec*DEG2RAD)*np.cos(hrang*DEG2RAD)
                        self.parang_f = np.arctan2(tmp1, tmp2)*RAD2DEG
        except KeyError:
            self.mjdtime = None
            self.lst_i = None
            self.az_i = None
            self.alt_i = None
            self.parang_i = None
            self.lst_f = None
            self.az_f = None
            self.alt_f = None
            self.parang_f = None
            self.error_id = 1
        try:
            self.tau_time = hdr['TAU_TIME']  # sec
        except KeyError:
            self.tau_time = None
            self.error_id = 1
        try:
            self.glon = hdr['GLON']  # deg
            self.glat = hdr['GLAT']
        except KeyError:
            if self.obs_ra is not None and self.obs_dec is not None:
                c = SkyCoord(self.obs_ra, self.obs_dec,
                             unit='deg', frame='fk5')
                self.glon = c.galactic.l.deg
                self.glat = c.galactic.b.deg
            else:
                self.glon = None
                self.glat = None
        try:
            self.square = hdr['SQUARE']
            if hdr['SQUARE'][-2:]==' A':
                self.square = hdr['SQUARE'][:-2] #slice off ' A'
        except KeyError:
            self.square = None
        try:
            self.az_star = hdr['AZ_STAR']
        except KeyError:
            self.az_star = None
        try:
            self.el_star = hdr['EL_STAR']
        except KeyError:
            self.el_star = None
        try:
            self.pa_star = hdr['PA_STAR']
        except KeyError:
            self.pa_star = None
        try:
            self.az_end = hdr['AZ_END']
        except KeyError:
            self.az_end = None
        try:
            self.el_end = hdr['EL_END']
        except KeyError:
            self.el_end = None
        try:
            self.pa_end = hdr['PA_END']
        except KeyError:
            self.pa_end = None

    def set_tsky(self, nside, skymap):
        """Sets the tsky attribute of the Image object, the
        brightness temp of the 341 MHz GSM map

        Parameters
        ----------
        nside : int
            Healpy nside parameter of skymap
        skymap : float array
            GSM map in healpy format
        """
        if self.glon is None or self.glat is None:
            self.tsky = None
        else:
            phi = self.glon
            theta = 90 - self.glat
            if phi < 0:
                phi += 360
            idx = hp.pixelfunc.ang2pix(nside, theta*DEG2RAD, phi*DEG2RAD)
            self.tsky = skymap[idx]

    def set_radius(self, scale):
        """Sets the radius attribute of the Image object.

        Parameters
        ----------
        scale : float
            Fraction between 0 and 1 of the image radius to use.
            The full size of the image field-of-view is multiplied
            by this number.
        """
        try:
            naxis = float(self.imsize.strip(')').split(',')[1])
            r = (naxis / 2.) * scale
            self.radius = round((r * self.pixel_scale) / 3600., 2)  # in deg
        except AttributeError:
            try:
                hdu, hdr = self.read()
                r = (hdr['NAXIS2'] / 2.) * scale
                self.radius = round(r * hdr['CDELT2'], 2)  # in deg
            except KeyError:
                self.radius = None

    def set_cycle(self, alwaysass):
        """Sets image config, cycle, semester, and ass_flag"""
        self.ass_flag = False
        for config in res_dict.keys():
            for cycle in res_dict[config].keys():
                if self.mjdtime <= res_dict[config][cycle]['mjd'][1] and self.mjdtime > res_dict[config][cycle]['mjd'][0]:
                    self.config = config
                    self.semester = res_dict[config][cycle]['semester']
                    self.cycle = cycle
                    break

        #set for VCSS snapshots
        if self.vcss:
            if self.noise > 30:
                self.ass_flag = False
            else:
                self.ass_flag = True
                
        # if ass_flag True, check if image within view of a "too-many-artifacts" bad source
        if self.ass_flag:
            bad_sources = {'3C286': SkyCoord(202.784583, 30.509167, unit='deg'),
                           '3C48': SkyCoord(24.422083, 33.159722, unit='deg'),
                           '3C147': SkyCoord(85.650417, 49.851944, unit='deg'),
                           '3C138': SkyCoord(80.291250, 16.639444, unit='deg'),
                           '3C84': SkyCoord(49.950417, 41.511667, unit='deg')}
            image_center = SkyCoord(self.obs_ra, self.obs_dec, unit='deg')
            bad_sep = 3.0  # if image within this dist of bad_source turn off ass_flag
            for src, loc in bad_sources.items():
                ang_sep = image_center.separation(loc).degree
                if ang_sep < bad_sep:
                    self.ass_flag = False
                    dbclasses_logger.info('SET_CYCLE WARNING: {} is in the '
                                          'field-of-view'.format(src))
                    break
        # set ass_flag if 'always associated' option enabled
        if alwaysass:
            self.ass_flag = True
        # unset ass_flag if a 'PER_FIELD' image
        if 'PER_FIEL' in self.filename:
            self.ass_flag = False

    def set_beam_image(self, pbdic, nobeamimage=False):
        """Calculates the primary beam image, offset & smeared, for the image.
        
        Parameters
        ----------
        pbdic : dictionary
           Primary beam dictionary
        nobeamimage : boolean
           Set True to skip beam image calculation (just parangs, zenith angles, weights)
        """
        if self.vcss is True:
            self.bmimg = beam_tools.Calc_Beam_Image_VCSS(self, pbdic, nobeamimage)
        else:
            self.bmimg = beam_tools.Calc_Beam_Image(self, pbdic, nobeamimage)

    def image_qa(self, params):
        """Performs quality checks on the image pre-source finding
        to flag images with the following issues:

        1. missing necessary header keywords (self.error_id = 1)
        2. number of visibilities (nvis) < min_nvis
        3. noise*sqrt(int. time) <= 0 or > max_sensitivity
        4. beam semi-major/semi-minor axis > max_ellip (too elliptical)
        5. target is NCP or a planet
        6. a known bright radio source is in the field-of-view (Cas A, Cygnus A, Taurus A, Hercules A, Virgo A (M87), Perseus A (3C84), the Sun, the Moon, Jupiter, Galactic Center)
        10. number of CLEAN iterations (niter) < min_niter
        11. BMIN in pixels < bpix_min or > bpix_max
        12. missing primary calibrators
        13. missing CLEAN components

        Images that fail checks 1-5,10-11 are aborted and do not proceed
        to source finding. Images are flagged if they fail check 6,
        but do continue on.

        """
        dbclasses_logger.info('Performing preliminary image quality checks...')

        #print('image_qa: error_id = ',self.error_id)
        # Check if image was missing any necessary header keywords
        if self.error_id == 1:
            dbclasses_logger.info('IMAGE FAILED QA: missing necessary header '
                                  'keyword(s).')
            return

        # Check number of visibilities (NVIS)
        min_nvis = params['min nvis']
        if self.nvis is not None and self.nvis < min_nvis:
            dbclasses_logger.info('IMAGE FAILED QA: number of visibilities '
                                  '(NVIS) {} < allowed minimum of {}.'.
                                  format(self.nvis, min_nvis))
            self.error_id = 2
            return
        else:
            pass

        # Check image sensitivity metric (noise * sqrt(int. time))
        max_sensitivity = params['max sensitivity metric']
        sensitivity_metric = self.noise * np.sqrt(self.tau_time)
        if sensitivity_metric <= 0:
            dbclasses_logger.info('IMAGE FAILED QA: sensitivity metric '
                                  '(noise x sqrt(int. time)) {} <= 0.'.
                                  format(sensitivity_metric))
            self.error_id = 3
            return
        if sensitivity_metric > max_sensitivity:
            dbclasses_logger.info('IMAGE FAILED QA: sensitivity metric '
                                  '(noise x sqrt(int. time)) {} > allowed '
                                  'maximum of {}.'.format(
                                      sensitivity_metric, max_sensitivity))
            self.error_id = 3
            return
        else:
            pass

        # Check number of CLEAN interations (NITER)
        min_niter = params['min niter']
        if self.niter is not None and self.niter < min_niter:
            dbclasses_logger.info('IMAGE FAILED QA: number of iterations '
                                  '(NITER) {} < allowed minimum of {}.'.
                                  format(self.niter, min_niter))
            self.error_id = 10
            return
        else:
            pass

        # Check beam ellipticity
        max_ellip = params['max beam axis ratio']
        axis_ratio = self.bmaj / self.bmin
        if axis_ratio > max_ellip:
            dbclasses_logger.info('IMAGE FAILED QA: beam axis ratio {} > '
                                  'allowed maximum of {}.'.format(
                                      axis_ratio, max_ellip))
            self.error_id = 4
            return
        else:
            pass

        # VLITE can't do planets or the NCP
        if self.obj == 'NCP' or self.obj == 'ncp':
            dbclasses_logger.info('IMAGE FAILED QA: pointing is at NCP.')
            self.error_id = 5
            return
        else:
            pass
        # Some planet observations have obs_ra, obs_dec = 0, 0
        if self.obs_ra == 0. and self.obs_dec == 0.:
            dbclasses_logger.info('IMAGE FAILED QA: planet observation with '
                                  'RA, Dec = 0, 0')
            self.error_id = 5
            return
        else:
            pass

        # Check BMIN in pixels
        bpix_min = params['min bpix']
        bpix_max = params['max bpix']
        bpix = self.bmin/self.pixel_scale
        if bpix < bpix_min:
            dbclasses_logger.info('IMAGE FAILED QA: bmin in pixels '
                                  '{} < allowed minimum of {}.'.
                                  format(bpix, bpix_min))
            self.error_id = 11
            return
        elif bpix > bpix_max:
            dbclasses_logger.info('IMAGE FAILED QA: bmin in pixels '
                                  '{} > allowed maximum of {}.'.
                                  format(bpix, bpix_max))
            self.error_id = 11
            return
        else:
            pass

        # Check primary calibrators
        if self.pri_cals is None:
            if self.filename.endswith('IMSC.fits') : #mosaics don't have
                pass
            else:
                self.error_id = 12
        else:
            pass

        # Check CLEAN components
        if self.cc is None:
            if self.filename.endswith('IMSC.fits') : #mosaics don't have
                pass
            else:
                self.error_id = 13
        else:
            pass

        # Check angular separation from problem sources
        sun = ephem.Sun()
        moon = ephem.Moon()
        jup = ephem.Jupiter()
        # convert from MJD to UTC ISO format
        jdt = Time(self.mjdtime, format='mjd', scale='utc')
        sun.compute(jdt.iso)
        moon.compute(jdt.iso)
        jup.compute(jdt.iso)
        bad_sources = {'Sun': SkyCoord(str(sun.a_ra), str(sun.a_dec),
                                       unit=('hourangle', 'deg')),
                       'Moon': SkyCoord(str(moon.a_ra), str(moon.a_dec),
                                        unit=('hourangle', 'deg')),
                       'Jupiter': SkyCoord(str(jup.a_ra), str(jup.a_dec),
                                           unit=('hourangle', 'deg')),
                       'Cas A': SkyCoord(350.866250, 58.811667, unit='deg'),
                       'Cen A': SkyCoord(201.365000, -43.019167, unit='deg'),
                       'Cyg A': SkyCoord(299.867917, 40.733889, unit='deg'),
                       'Her A': SkyCoord(252.783750, 4.992500, unit='deg'),
                       'Orion A': SkyCoord(83.818750, -5.389722, unit='deg'),
                       'Per A': SkyCoord(49.950417, 41.511667, unit='deg'),
                       'Tau A': SkyCoord(83.633333, 22.014444, unit='deg'),
                       'Virgo A': SkyCoord(187.705833, 12.391111, unit='deg'),
                       'GC': SkyCoord(266.416833, -29.007806, unit='deg')}

        image_center = SkyCoord(self.obs_ra, self.obs_dec, unit='deg')

        min_sep = 999999.99
        for src, loc in bad_sources.items():
            ang_sep = image_center.separation(loc).degree
            if src == 'Sun':
                self.sunsep = ang_sep
            while ang_sep < min_sep:
                min_sep = ang_sep
                self.nearest_problem = src
                self.separation = ang_sep
        if min_sep <= self.radius:
            dbclasses_logger.info('IMAGE QA WARNING: {} is in the '
                                  'field-of-view'.format(self.nearest_problem))
            self.error_id = 6
            return
        else:
            pass

        # all checks passed:
        dbclasses_logger.info('...image passed.')

    def source_qa(self, sources, params):
        """Approximates the expected number of sources in an
        image based on source counts from WENSS scaled by
        noise and fitted beam. See EP's PythonTools.

        """
        max_src_metric = params['max source count metric']

        dbclasses_logger.info('Performing source count quality checks...')
        # First, check if there are any sources
        if len(sources) < 1:
            dbclasses_logger.info('IMAGE FAILED QA: no sources were detected.')
            self.error_id = 8
            return
        else:
            pass

        # Get number of actual sources within central 1.5 degrees
        nsrc_cut = 0
        for src in sources:
            if src.dist_from_center <= 1.5:
                nsrc_cut += 1

        '''
        # Estimate number of expected sources within central 1.5 degrees
        nsrc_exp = beam_tools.expected_nsrc(self.pri_freq, self.noise)
        # Compute metric for determining if source count is way off
        nsrc_metric = (float(nsrc_cut) - nsrc_exp) / nsrc_exp

        if nsrc_metric > max_src_metric:
            dbclasses_logger.info('IMAGE FAILED QA: source count metric '
                                  '{} > allowed max {}.'.format(
                                      nsrc_metric, max_src_metric))
            self.error_id = 9
            return
        else:
            pass
        '''

        dbclasses_logger.info('...image passed.')

    def log_attrs(self, pybdsfdir):
        """Image object method to assign values to object 
        attributes which come from source finding using info 
        in the PyBDSF output log file. This is useful when 
        source finding has already been run on an image and the 
        results are being recorded after the fact using the 
        PyBDSF output files.

        Parameters
        ----------
        pybdsfdir : str
            Directory path to the PyBDSF output files.

        Returns
        -------
        rms_box : str
            PyBDSF parameter used during source finding.
        gresid_std : float
            Standard deviation of the background in the residual 
            image after removing gaussian fitted sources.
        raw_rms : float
            Estimated noise in the image before source extraction.
        """
        prefix = re.findall('.*\/(.*)', img.filename)[0]
        logname = prefix + '.pybdsf.log'
        log = os.path.join(pybdsfdir, logname)
        with open(path) as f:
            log = f.read()
        x = re.findall('.*rms_box.*:\s(\(.*\))', log)
        rms_box = x[-1]
        y = re.findall('.*std. dev:\s(.*)\s\(', log)
        gresid_std = float(y[-1]) * 1000.  # mJy/beam
        z = re.findall('.*raw rms =\s+(.*) mJy', log)
        raw_rms = float(z[-1])
        return rms_box, gresid_std, raw_rms


class DetectedSource(object):
    """Class of objects to store elliptical Gaussian fit 
    properties of sources found and measured by PyBDSF. 
    Attribute values are translated from the class of objects 
    output by PyBDSF.

    Attributes
    ----------
    src_id : int
        Uniquely identifies the source in a given image.
    isl_id : int
        Uniquely identifies an island in a given image from
        which sources are formed.
    image_id : int
        Uniquely identifies the image the source comes from.
    ra : float
        Right ascension (degrees).
    e_ra : float
        Error on the right ascension (degrees).
    dec : float
        Declination (degrees).
    e_dec : float
        Error on the declination (degrees).
    total_flux : float
        Total integrated flux (mJy).
    e_total_flux : float
        Error on the total integrated flux (mJy).
    peak_flux : float
        Peak flux density per beam (mJy/beam).
    e_peak_flux : float
        Error on the peak flux (mJy/beam).
    ra_max : float
        Right ascension of the source maximum brightness (degrees).
    e_ra_max : float
        Error on the right ascension of the maximum (degrees).
    dec_max : float
        Declination of the source maximum brightness (degrees).
    e_dec_max : float
        Error on the declination of the maximum (degrees).
    maj : float
        FWHM of the source major axis (arcsec).
    e_maj : float
        Error on the source major axis size (arcsec).
    min : float
        FWHM of the source minor axis (arcsec).
    e_min : float
        Error on the source minor axis size (arcsec).
    pa : float
        Position angle of the source major axis measured east
        of north (degrees).
    e_pa : float
        Error on the source position angle (degrees).
    dc_maj : float
        FWHM of the deconvolved major axis (arcsec).
    e_dc_maj : float
        Error on the deconvolved major axis (arcsec).
    dc_min : float
        FWHM of the deconvolved minor axis (arcsec).
    e_dc_min : float
        Error on the deconvolved minor axis (arcsec).
    dc_pa : float
        Position angle of the deconvolved major axis measured
        east of north (degrees).
    e_dc_pa : float
        Error on the deconvolved position angle (degrees).
    total_flux_isl : float
        Total integrated flux density of the island (mJy).
    total_flux_islE : float
        Error on the total integrated flux density of the island (mJy).
    rms_isl : float
        Average background rms of the island (mJy/beam).
    mean_isl : float
        Average background mean of the island (mJy/beam).
    resid_rms : float
        Average residual background rms of the island (mJy/beam).
    resid_mean : float
        Average residual background mean of the island (mJy/beam).
    code : str
        Defines source structure: 'S' = isolated single-Gaussian source,
        'C' = single-Gaussian source in island with other sources,
        'M' = multi-Gaussian source.
    assoc_id : int
        Links to associated source in the database **assoc_source** table.
    dist_from_center : float
        Angular separation between the source location and image
        pointing center (degrees).
    polar_angle : float
        Polar angle (East of North) of source in image (degrees).
    compactness: float
        Measure of source extent. 
        >= 1 : point-like
        < 1  : extended
    clean : boolean
        True if source was CLEANed.
    xpix : float
        X pixel value of source
    ypix : float
        Y pixel value of source
    id : int
        Uniquely identifies the source in the **assoc_source** table.
    res_class : str
        Resolution class of the image in which the source was found.
    ndetect : int
        Number of times this same source has been detected in other
        VLITE images.
    ns : int
        Number of 'S' detections of this source
    nc : int
        Number of 'C' detections of this source
    nm : int
        Number of 'M' detections of this source
    nmatches : int
        Number of catalogs which contain this source.
    detected : bool
        Whether or not the source was detected in the image (used
        for the **vlite_unique** table).
    ave_total : float
        Weighted average of corrected total flux for associated sources 
        observed more than once (mJy)
    e_ave_total : float
        Error on the weighted average of corrected total flux (mJy)
    ave_peak : float
        Weighted average of corrected peak flux for associated sources 
        observed more than once (mJy/beam)
    e_ave_peak : float
        Error on the weighted average of corrected peak flux (mJy/beam)
    v_total : float
        Variability metric of associated source ave_total
    v_peak : float
        Variability metric of associated source ave_peak
    eta_total : float
        Variability significance metric of associated source ave_total
    eta_peak : float
        Variability significance metric of associated source ave_peak
    nn_src_id : int
        src_id of nearest neighbor in image
    nn_dist : float
        Distance to nearest neighbor (arcsec)

    References
    ----------
    Refer to the PyBDSF documentation[1]_ for definitions 
    of their output columns.

    .. [1] http://www.astron.nl/citt/pybdsm/write_catalog.html#definition-of-output-columns
    """

    def __init__(self):
        self.src_id = None
        self.isl_id = None
        self.image_id = None
        self.ra = None
        self.e_ra = None
        self.dec = None
        self.e_dec = None
        self.total_flux = None
        self.e_total_flux = None
        self.peak_flux = None
        self.e_peak_flux = None
        self.ra_max = None
        self.e_ra_max = None
        self.dec_max = None
        self.e_dec_max = None
        self.maj = None
        self.e_maj = None
        self.min = None
        self.e_min = None
        self.pa = None
        self.e_pa = None
        self.dc_maj = None
        self.e_dc_maj = None
        self.dc_min = None
        self.e_dc_min = None
        self.dc_pa = None
        self.e_dc_pa = None
        self.total_flux_isl = None
        self.total_flux_islE = None
        self.rms_isl = None
        self.mean_isl = None
        self.resid_rms = None
        self.resid_mean = None
        self.code = None
        self.assoc_id = None
        self.dist_from_center = None
        self.polar_angle = None
        self.compactness = None
        self.clean = None
        self.xpix = None
        self.ypix = None
        self.nn_src_id = None
        self.nn_dist = None
        # assoc_source attributes
        self.id = None
        self.res_class = None
        self.ndetect = None
        self.ns = None
        self.nc = None
        self.nm = None
        self.nmatches = None
        self.ave_total = None
        self.e_ave_total = None
        self.ave_peak = None
        self.e_ave_peak = None
        self.v_total = None
        self.v_peak = None
        self.eta_total = None
        self.eta_peak = None
        # vlite_unique attribute
        self.detected = None

    def cast(self, origsrc):
        """Attributes of a ``bdsf.gaul2srl.Source`` object
        output from ``bdsf.process_image`` are re-cast to define
        attributes of a DetectedSource object.

        Parameters
        ----------
        origsrc : ``bdsf.gaul2srl.Source`` instance
            Object which stores all the measured properties
            of a source extracted from an image in PyBDSF.

        """
        self.src_id = origsrc.source_id
        self.isl_id = origsrc.island_id
        self.ra = origsrc.posn_sky_centroid[0]  # deg
        self.e_ra = origsrc.posn_sky_centroidE[0]  # deg
        self.dec = origsrc.posn_sky_centroid[1]  # deg
        self.e_dec = origsrc.posn_sky_centroidE[1]  # deg
        self.total_flux = origsrc.total_flux * 1000.  # mJy
        self.e_total_flux = origsrc.total_fluxE * 1000.  # mJy
        self.peak_flux = origsrc.peak_flux_centroid * 1000.  # mJy/beam
        self.e_peak_flux = origsrc.peak_flux_centroidE * 1000.  # mJy/beam
        self.ra_max = origsrc.posn_sky_max[0]  # deg
        self.e_ra_max = origsrc.posn_sky_maxE[0]  # deg
        self.dec_max = origsrc.posn_sky_max[1]  # deg
        self.e_dec_max = origsrc.posn_sky_maxE[1]  # deg
        self.maj = origsrc.size_sky[0] * 3600.  # arcsec
        self.e_maj = origsrc.size_skyE[0] * 3600.  # arcsec
        self.min = origsrc.size_sky[1] * 3600.  # arcsec
        self.e_min = origsrc.size_skyE[1] * 3600.  # arcsec
        self.pa = origsrc.size_sky[2]  # deg
        self.e_pa = origsrc.size_skyE[2]  # deg
        self.dc_maj = origsrc.deconv_size_sky[0] * 3600.  # deconvolved
        self.e_dc_maj = origsrc.deconv_size_skyE[0] * 3600.
        self.dc_min = origsrc.deconv_size_sky[1] * 3600.
        self.e_dc_min = origsrc.deconv_size_skyE[1] * 3600.
        self.dc_pa = origsrc.deconv_size_sky[2]
        self.e_dc_pa = origsrc.deconv_size_skyE[2]
        self.total_flux_isl = origsrc.total_flux_isl * 1000.  # mJy
        self.total_flux_islE = origsrc.total_flux_islE * 1000.  # mJy
        self.rms_isl = origsrc.rms_isl * 1000.  # mJy/beam
        self.mean_isl = origsrc.mean_isl * 1000.  # mJy/beam
        self.resid_rms = origsrc.gresid_rms * 1000.  # mJy/beam
        self.resid_mean = origsrc.gresid_mean * 1000.  # mJy/beam
        self.code = origsrc.code

    def calc_center_dist(self, imobj):
        """Calculates a detected source's angular distance
        from the image center in degrees and sets
        the ``dist_from_center`` attribute.

        Parameters
        ----------
        imobj : ``database.dbclasses.Image`` instance
            Initialized Image object with attribute values
            set from header info.
        """
        img_center = SkyCoord(imobj.obs_ra, imobj.obs_dec, unit='deg')
        src_loc = SkyCoord(self.ra, self.dec, unit='deg')
        self.dist_from_center = img_center.separation(src_loc).degree

    def calc_polar_angle(self, imobj):
        """Calculates a source's polar angle, E of N
        and sets the ``polar_angle`` attribute.

        Parameters
        ----------
        imobj : ``database.dbclasses.Image`` instance
            Initialized Image object with attribute values
            set from header info.
        """
        theta = np.arctan2(imobj.xref-self.xpix, self.ypix-imobj.yref)  # rad
        self.polar_angle = theta * 180.0/np.pi  # deg

    def calc_pixel_coords(self, imobj):
        """Calculates a detected source's position in pixels
        and sets the ``xpix`` & ``ypix`` attributes.

        Parameters
        ----------
        imobj : ``database.dbclasses.Image`` instance
            Initialized Image object with attribute values
            set from header info.
        """
        self.xpix, self.ypix = beam_tools.EQtoPIX(self.ra, self.dec, imobj.wcsobj)


    def correct_flux(self, imobj, pbdic):
        """Uses primary beam image to correct all flux measurements
        from PyBDSF.

        Parameters
        ----------
        imobj : ``database.dbclasses.Image`` instance
            Object with image attributes
        """
        # Correct for beam response
        # Use source pixel coords to get pixel value in beam image
        pb_power = imobj.bmimg[int(self.ypix),int(self.xpix)]
        pb_err = pbdic[imobj.pbkey].error
        
        # List of attributes to correct
        attrs = ['total_flux', 'peak_flux', 'total_flux_isl', 'rms_isl',
                 'mean_isl', 'resid_rms', 'resid_mean', 'e_total_flux',
                 'e_peak_flux', 'total_flux_islE']
        for attr in attrs:
            corr_val = getattr(self, attr) / pb_power
            setattr(self, attr, corr_val)
        # Add systematic uncertainty to all flux errors
        #self.e_total_flux = np.sqrt(
        #    (self.total_flux * pb_err)**2. + self.e_total_flux**2.)
        #self.e_peak_flux = np.sqrt(
        #    (self.peak_flux * pb_err)**2. + self.e_peak_flux**2.)
        #self.total_flux_islE = np.sqrt(
        #    (self.total_flux_isl * pb_err)**2. + self.total_flux_islE**2.)

    def calc_snr(self):
        # Compute S/N of source detection
        self.snr = self.peak_flux / self.rms_isl

    # nulls need special handling
    def correct_flux_null(self, pri_freq):
        """Applies a 1-D radial correction to null flux measurements
        from PyBDSF.

        Parameters
        ----------
        pri_freq : float
            Primary frequency of the observations in GHz.
        """
        # Correct for beam response - only 1D (sym. beam) for now
        pb_power, pb_err = beam_tools.find_nearest_pbcorr(self.dist_from_center,
                                                          pri_freq)
        # List of attributes to correct
        attrs = ['total_flux', 'peak_flux', 'e_total_flux', 'e_peak_flux']
        for attr in attrs:
            corr_val = getattr(self, attr) / pb_power
            setattr(self, attr, corr_val)
        # Compute S/N of source detection
        #self.snr = (self.peak_flux - self.mean_isl) / self.rms_isl
        # Add systematic uncertainty to all flux errors
        #self.e_total_flux = np.sqrt(
        #    (self.total_flux * pb_err)**2. + self.e_total_flux**2.)
        #self.e_peak_flux = np.sqrt(
        #    (self.peak_flux * pb_err)**2. + self.e_peak_flux**2.)
        # self.total_flux_islE = np.sqrt(
        #    (self.total_flux_isl * pb_err)**2. + self.total_flux_islE**2.)

    def calc_compactness(self, imobj):
        """Calculates a detected source's compactness
        from its flux ratio, SNR, and array config.
        Sets the ``compactness`` attribute.

        Parameters
        ----------
        imobj : ``database.dbclasses.Image`` instance
            Initialized Image object with attribute values
            set from header info.
        """
        self.compactness = None
        flag = False
        if imobj.mosaic:
            pass
        elif imobj.vcss:
            config = 'VCSS'
            if imobj.mjdtime > 58000 and imobj.mjdtime < 58200:
                cycle = '1.1'
                flag = True
            elif imobj.mjdtime > 58500 and imobj.mjdtime < 58700:
                cycle = '1.2'
                flag = True
            elif imobj.mjdtime > 59040 and imobj.mjdtime < 59160:
                cycle = '2.1'
                flag = True
            elif imobj.mjdtime > 59400 and imobj.mjdtime < 59700:
                cycle = '2.2'
                flag = True
            elif imobj.mjdtime > 59900 and imobj.mjdtime < 60150:
                cycle = '3.1'
                flag = True
            else:
                pass
        else:
            for cfg in compactness_dict.keys():
                if cfg == imobj.config:
                    for cyc in compactness_dict[imobj.config].keys():
                        if cyc == imobj.cycle:
                            config = imobj.config
                            cycle = imobj.cycle
                            flag = True
                            break
        if flag:
            a0 = compactness_dict[config][cycle]['a0']
            a1 = compactness_dict[config][cycle]['a1']
            c0 = compactness_dict[config][cycle]['c0']
            c1 = compactness_dict[config][cycle]['c1']
            cfit = a0+sqrt(a1 + c0*pow(self.snr, c1))
            ratio = self.total_flux/self.peak_flux
            self.compactness = cfit/ratio


def dict2attr(obj, dictionary):
    """Converts dictionary key-value pairs to object attributes.

    Parameters
    ----------
    obj : class instance
        Any object initialized from a class.
    dictionary : dict
        Dictionary key-value pairs become attribute-value pairs.
    """
    for key in dictionary.keys():
        setattr(obj, key, dictionary[key])


def translate_from_txtfile(img, pybdsfdir):
    """Creates a list of DetectedSource objects
    from reading in a source list text file output
    by PyBDSF.

    Parameters
    ----------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute
        values set from header info.
    pybdsfdir : str
        Directory path to PyBDSF files.

    Returns
    -------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute
        values updated with source finding results.
    sources : list of DetectedSource objects
        DetectedSource objects with attribute values
        set from PyBDSF output source list text file.
    """
    # Set the file path
    prefix = re.findall('.*\/(.*)\.', img.filename)[0]
    srl = prefix + '.pybdsm.srl'
    catalog = os.path.join(pybdsfdir, srl)
    # Read the catalog
    dbclasses_logger.info('Extracting sources from {}.'.format(srl))
    sources = pybdsfcat.read_catalog(catalog)

    # Count the number of sources
    img.nsrc = len(sources)

    # Extract rms_box size and noise from log file
    img.rms_box, img.noise, raw_rms = img.log_attrs(pybdsfdir)

    return img, sources


def translate(img, out):
    """Method to translate PyBDSF output within the
    pipeline to DetectedSource objects and update 
    the Image object.

    Parameters
    ----------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info.
    out : ``bdsf.image.Image`` instance
        The object output by PyBDSF after running its
        source finding task ``process_image()``. Contains
        a list of ``bdsf.gaul2srl.Source`` objects which
        are translated into DetectedSource objects.

    Returns
    -------
    newsrcs : list
        List of ``database.dbclasses.DetectedSource`` objects.
        Attributes of each object are set from the PyBDSF
        output object.    
    """
    # Translate PyBDSF output source objects to DetectedSource objects
    newsrcs = []
    print('**************************************')
    print('len sources = ',len(out.sources))
    for oldsrc in out.sources: #flagged sources are in out.dsources (d for dummy)
        newsrcs.append(DetectedSource())
        newsrcs[-1].cast(oldsrc)
        newsrcs[-1].image_id = img.id
        newsrcs[-1].calc_center_dist(img)
        newsrcs[-1].calc_pixel_coords(img)
        newsrcs[-1].calc_polar_angle(img)
        #initialize to 0
        newsrcs[-1].ns = 0
        newsrcs[-1].nc = 0
        newsrcs[-1].nm = 0
        if oldsrc.code == 'S': newsrcs[-1].ns = 1
        elif oldsrc.code == 'C': newsrcs[-1].nc = 1
        elif oldsrc.code == 'M': newsrcs[-1].nm = 1
        else: print('ERROR! CODE'+oldsrc.code+'NOT RECOGNIZED!')
    '''
    ############################## works, but not being written to database!:
    print('len dsources =',len(out.dsources))
    for oldsrc in out.dsources: #flagged sources are in out.dsources (d for dummy)
        newsrcs.append(DetectedSource())
        newsrcs[-1].cast(oldsrc)
        newsrcs[-1].image_id = img.id
        newsrcs[-1].calc_center_dist(img)
        newsrcs[-1].calc_pixel_coords(img)
        newsrcs[-1].calc_polar_angle(img)
        #initialize to 0
        newsrcs[-1].ns = 0
        newsrcs[-1].nc = 0
        newsrcs[-1].nm = 0
        if oldsrc.code == 'S': newsrcs[-1].ns = 1
        elif oldsrc.code == 'C': newsrcs[-1].nc = 1
        elif oldsrc.code == 'M': newsrcs[-1].nm = 1
        else: print('ERROR! CODE'+oldsrc.code+'NOT RECOGNIZED!')
    print('len newsources = ',len(newsrcs))
    print('*****************************************')
    ##############################
    '''
    return newsrcs


# nulls just use pybdsf islands and require
# different handling than sources
def translate_null(img, out, coords):
    """Method to translate PyBDSF islands output within the
    pipeline to DetectedSource objects 

    Parameters
    ----------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info.
    out : ``bdsf.image.Image`` instance
        The object output by PyBDSF after running its
        source finding task ``process_image()``. Contains
        a list of ``bdsf.gaul2srl.Island`` objects which
        are translated into DetectedSource objects.
    coords : list of (RA, Dec) tuples specfying Island coords

    Returns
    -------
    newsrcs : list
        List of ``database.dbclasses.DetectedSource`` objects.
        Attributes of each object are set from the PyBDSF
        output object.    
    """
    # Translate PyBDSF output island objects to DetectedSource objects
    newsrcs = []
    for n, oldsrc in enumerate(out.islands):
        newsrcs.append(DetectedSource())
        newsrcs[-1].ra = coords[n][0]
        newsrcs[-1].dec = coords[n][1]
        newsrcs[-1].image_id = img.id
        # set total flux to island max value,
        # this gave best results compared to PySE force-fitting
        newsrcs[-1].total_flux = oldsrc.max_value*1000.
        # set flux err to island total flux err
        newsrcs[-1].e_total_flux = oldsrc.total_fluxE*1000.
        # set peak to total
        newsrcs[-1].peak_flux = oldsrc.max_value*1000.
        newsrcs[-1].e_peak_flux = oldsrc.total_fluxE*1000.
        newsrcs[-1].calc_center_dist(img)
        newsrcs[-1].calc_polar_angle(img)

    return newsrcs


def set_fromnx(img, smear_time):
    """Counts number of SN tables, number of observing 
    intervals; set list of times and weigths for primary 
    beam calculation from NX table in UVOUT file    

    Parameters
    ----------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute
        values set from header info.
    smear_time : Max time step for smearing beam [s]

    Returns
    -------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute
        nsn updated.
    """
    # Skip for VCSS images
    # (don't need times or weights for beam calc)
    if img.vcss:
        img.smeartime = 2.0
        img.pbparangs = []
        img.pbweights = []
        img.pbza = []
        img.pbtimes = []
        img.nsn = None
        img.pb_flag = True
        return img

    # Skip for VCSS mosaics
    # (don't need times or weights for beam calc)
    if img.mosaic:
        img.smeartime = None
        img.pbparangs = []
        img.pbweights = []
        img.pbza = []
        img.pbtimes = []
        img.nsn = None
        img.pb_flag = False
        return img

    img.smeartime = smear_time
    
    # Determine dir with uvout file
    a = img.filename.split('/')
    uvdir = ''
    for i in range(len(a)-2):
        uvdir += (a[i]+'/')
    uvdir += 'UVFiles/'
    # Check dir
    if not os.path.isdir(uvdir):
        img.nsn = None
        img.pb_flag = False
        dbclasses_logger.info('Did not find UVFiles dir {}.'.format(uvdir))
        dbclasses_logger.info('Cannot calc primary beam for this image')
        return img

    # Determine name of uvout file
    b = a[-1].split('.')
    uvname1 = ''
    for i in range(len(b)-2):
        uvname1 += (b[i]+'.')
    # for time chop images:
    if img.filename.endswith('ITime.fits'):
        uvname = uvname1 + 'UVOUT.fits'
    # for dailies:
    else:
        uvname = uvname1 + 'uvout.gz'

    # Check file
    if not os.path.isfile(os.path.join(uvdir, uvname)):
        # unzipped?
        uvname = uvname1 + 'uvout'
        if not os.path.isfile(os.path.join(uvdir, uvname)):
            img.nsn = None
            img.pb_flag = False
            dbclasses_logger.info(
                'Did not find uvout file, with or w/o .gz {}'.format(uvname))
            dbclasses_logger.info(
                'Cannot calc primary beam for this image')
            return img

    # Count number of SN tables & set times, weights from NX table
    img.nsn = 0
    img.nbeam = 0
    img.pb_flag = True
    with fits.open(os.path.join(uvdir, uvname)) as hdu:
        for i in range(1, len(hdu)):
            if hdu[i].header['EXTNAME'] == 'AIPS SN':
                img.nsn += 1
            if hdu[i].header['EXTNAME'] == 'AIPS NX':
                img.pbweights = []
                img.pbtimes = []
                img.pbparangs = []
                img.pbza = []
                hdr = hdu[i].header
                img.ninterval = hdr['NAXIS2']
                data = hdu[i].data
                img.max_dt = 86400*np.max(data['TIME INTERVAL']) #s
                #'TIME' is center time of interval in days since reference time
                ti = data['TIME']-(data['TIME INTERVAL']/2)-img.startime-(1./86400) #day
                tf = data['TIME']+(data['TIME INTERVAL']/2)-img.startime+(1./86400) #day
                dt = tf-ti #day
                dvis = data['END VIS']-data['START VIS']+1 #includes end pts
                img.nvisnx = int(np.sum(dvis))
                #sampling includes end pts
                for j in range(img.ninterval):
                    npb = ceil(dt[j]/(smear_time/86400)) + 1
                    dsmeartime = dt[j]/(npb-1)
                    for k in range(npb):
                        img.pbweights.append((dvis[j]/img.nvisnx)/npb)
                        img.pbtimes.append(ti[j]+(k*dsmeartime)) #day
                        img.nbeam += 1
                '''
                #sampling starts at half step - including end points is better
                for j in range(img.ninterval):
                    npb = ceil(dt[j]/(smear_time/86400))
                    dsmeartime = dt[j]/npb
                    for k in range(npb):
                        img.pbweights.append((dvis[j]/img.nvisnx)/npb)
                        img.pbtimes.append(ti[j]+((k+0.5)*dsmeartime)) #day
                        img.nbeam += 1
                '''
    if img.nbeam == 0: #No NX table found
        img.nbeam = None
        img.pb_flag = False
        img.ass_flag = False
        if img.error_id == None: #don't overwrite previously set error
            img.error_id = 14 
    print(uvname, img.nsn, img.nbeam, img.ninterval, img.max_dt, img.nvisnx, np.sum(img.pbweights))
    dbclasses_logger.info('{}: nsn= {}; nbeam= {}; ninterval= {}; max_dt= {}; nvisnx= {}; sumweights= {}'.format(uvname, img.nsn, img.nbeam, img.ninterval, img.max_dt, img.nvisnx, np.sum(img.pbweights)))
    
    return img


def init_image(impath, alwaysass, smear_time):
    """Initializes an object of the Image class and sets values 
    for its attributes from the fits file header using
    the ``header_attrs`` object method.

    Parameters
    ----------
    alwaysass : 'always associate' option. If True image sources
              will always be associated
    smear_time : 'smear time' option. Determines number of beam
              samples in each observing block
    """
    imobj = Image()
    imobj.process_image(impath)
    hdu, hdr = imobj.read()
    # Use header info to set attributes
    imobj.header_attrs(hdr)
    # Fix header keywords for PyBDSF, if necessary
    header_changed = False
    if hdr['NAXIS'] > 2:
        if hdr['CTYPE3'] == 'SPECLNMF':
            hdr['CTYPE3'] = 'FREQ'
            header_changed = True
    try:
        hdr['BMAJ']
    except KeyError:
        hdr['BMAJ'] = imobj.bmaj / 3600.  # deg
        hdr['BMIN'] = imobj.bmin / 3600.
        hdr['BPA'] = imobj.bpa
        header_changed = True
    if header_changed:
        imobj.write(hdu, hdr)

    # Set cycle, semester, ass_flag
    imobj.set_cycle(alwaysass)

    # Set attributes from NX table of uvout
    imobj = set_fromnx(imobj, smear_time)

    # Read primary calibrators from history extension
    pri_cals = []
    for i in range(1, len(hdu)):
        if hdu[i].header['EXTNAME'] == 'History':
            for j in hdu[i].data['ENTRY']:
                if 'Primary Calibrator' in j:
                    k = j.split("=")[1].split("M")[0].lstrip().rstrip()
                    if k not in pri_cals:
                        pri_cals.append(k)
    if len(pri_cals) > 0:
        imobj.pri_cals = pri_cals
    else:
        imobj.pri_cals = None

    # Read CLEAN components from CC extension
    xo = hdr['CRPIX1']
    yo = hdr['CRPIX2']
    w = wcs.WCS(hdr)
    cc = []
    for i in range(1, len(hdu)):
        if 'CC' in hdu[i].header['EXTNAME']:
            for j in range(len(hdu[i].data['DELTAX'])):
                x = xo+(hdu[i].data['DELTAX'][j]/hdr['CDELT1'])
                y = yo+(hdu[i].data['DELTAY'][j]/hdr['CDELT2'])
                pixcrd = w.wcs_pix2world([[x, y, 1, 1]], 1)
                ra, dec = pixcrd[0][0], pixcrd[0][1]
                cc.append((ra, dec))
    if len(cc) > 0:
        imobj.cc = cc
    else:
        imobj.cc = None

    hdu.close()
    return imobj



def getimgmjd(impath):
    """Initializes an object of the ImgMjd class and sets values 
    for its attributes from the fits file header using
    the ``init_image`` object method.
    """
    imobj = ImgMjd()
    imobj.init_image(impath)
    return imobj
