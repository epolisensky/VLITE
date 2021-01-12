"""This module defines classes and their methods for
Image and DetectedSource objects.

"""
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
from math import sqrt
from datetime import datetime

########################
# temporarily needed while USNO sites are down for modernization
#  expected completeion 30 Apr 2020
from astropy.utils import iers
iers.conf.iers_auto_url = u'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
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
#  the acceptable range of image bmins for source association
res_dict = {
    'A': {'1': {'mjd': [58179, 58280], 'bmin': [2.85, 4.27], 'semester': '2018A'},
          '2': {'mjd': [58697, 58778], 'bmin': [3.26, 4.90], 'semester': '2019A'}},
    'BnA': {'1': {'mjd': [58148, 58179], 'bmin': [0, 0], 'semester': ''},
            '2': {'mjd': [58660, 58697], 'bmin': [0, 0], 'semester': ''}},
    'B': {'1': {'mjd': [57996, 58148], 'bmin': [10.77, 16.16], 'semester': '2017B'},
          '2': {'mjd': [58534.24, 58660], 'bmin': [10.18, 15.28], 'semester': '2019A'},
          '3': {'mjd': [59026, 99999], 'bmin': [11.65, 17.47], 'semester': '2020A'}},
    'CnB': {'1': {'mjd': [57994, 57996], 'bmin': [0, 0], 'semester': ''},
            '2': {'mjd': [58519, 58534.24], 'bmin': [0, 0], 'semester': ''},
            '3': {'mjd': [59008, 59026], 'bmin': [0, 0], 'semester': ''}},
    'C': {'V10': {'mjd': [56986, 57953.99], 'bmin': [32.2, 65.91], 'semester': 'many'},
          '1': {'mjd': [57954, 57994], 'bmin': [43.94, 65.91], 'semester': '2017A'},
          '2': {'mjd': [58441, 58519], 'bmin': [32.20, 48.30], 'semester': '2018B'},
          '3': {'mjd': [58885, 59008], 'bmin': [39.69, 59.54], 'semester': '2020A'}},
    'DnC': {'3': {'mjd': [58876, 58885], 'bmin': [0, 0], 'semester': ''}},
    'D': {'2': {'mjd': [58360, 58440], 'bmin': [117.58, 176.37], 'semester': '2018A'},
          '3': {'mjd': [58802, 58876], 'bmin': [133.07, 199.60], 'semester': '2019B'}},
    'AnB': {'V10': {'mjd': [56986, 57953.99], 'bmin': [0, 0], 'semester': 'many'}}
}

# Define beam offset & dictionary for each primary observing frequency.
# Beams are offset from phase center by amount "offset" in given direction.
# Angles are in deg, "E" of zenith direction (counter-clkwise looking into dish)
# *Note P band primary observing images do NOT need offset correction
offset = 6.5/60  # deg
offset_dict = {'1.5': -174.1,
               '3': 11.6,
               '6': 75.2,
               '10': 113.7,
               '15': -42.4,
               '22': -64.1,
               '33': -106.9,
               '45': -85.5
               }


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
        Right ascension of image pointing center (degrees).
    obs_dec : float
        Declination of image pointing center (degrees).
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
    az_f : float
        Azimuth of image pointing center at mjdtime+duration from astropy (degrees).
    alt_f : float
        Altitude of image pointing center at mjdtime+duration from astropy (degrees).
    parang_f : float
        Parallactic angle of image pointing center at mjdtime+duration from astropy (degrees).
    lst : float
        Local Sidereal Time of image at start time (hours) 
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
        Integration time on source in seconds.
    duration : float
        Total time spanned by the observations in seconds.
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
        self.lst = None
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
        # start with priband
        try:
            priband = hdr['PRIBAND']
            if float(priband[:-3]) > 100:
                self.pri_freq = float(priband[:-3])/1000  # GHz
            else:
                self.pri_freq = float(priband[:-3])  # GHz
        except:
            # VCSS mosaics/snapshots
            if self.filename.endswith('IMSC.fits') or 'VCSS' in self.filename:
                self.pri_freq = 3
            else:
                pri = re.findall('\/([0-9.]+[A-Z])', self.filename)
                try:
                    if pri[0][-1:] == 'M':
                        self.pri_freq = float(pri[0][:-1]) / 1000.  # GHz
                    else:
                        self.pri_freq = float(pri[0][:-1])  # GHz
                    # project code 13B-266 is at 1.5 GHz
                    if self.pri_freq == 13:
                        self.pri_freq = 1.5
                except IndexError:
                    self.pri_freq = None
        ###
        try:
            naxis1 = hdr['NAXIS1']
            naxis2 = hdr['NAXIS2']
            self.imsize = str((naxis1, naxis2))  # pixels
        except KeyError:
            self.imsize = None
            self.error_id = 1

        self.wcsobj = wcs.WCS(hdr).celestial

        try:
            self.obs_ra = hdr['OBSRA']  # deg
        except KeyError:
            self.obs_ra = None
            self.error_id = 1
        if self.obs_ra == 0 or self.obs_ra is None:
            try:
                self.obs_ra = hdr['CRVAL1']  # deg
                self.error_id = None
            except KeyError:
                self.obs_ra = None
                self.error_id = 1
        try:
            self.obs_dec = hdr['OBSDEC']  # deg
        except KeyError:
            self.obs_dec = None
            self.error_id = 1
        if self.obs_dec == 0 or self.obs_dec is None:
            try:
                self.obs_dec = hdr['CRVAL2']  # deg
                self.error_id = None
            except KeyError:
                self.obs_dec = None
                self.error_id = 1
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
            self.config = hdr['CONFIG']
        except KeyError:
            self.config = None
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
                except:
                    self.mjdtime = Time(datetime(int(date[0]), int(date[1]), int(
                        date[2]))).mjd  # day
                t = Time(self.mjdtime, format='mjd')
                self.lst = t.sidereal_time('apparent', VLALON).hour  # hrs
                if self.obs_ra is not None and self.obs_dec is not None:
                    coord = SkyCoord(self.obs_ra, self.obs_dec,
                                     unit='deg', frame='fk5')
                    altaz = coord.transform_to(
                        AltAz(obstime=t, location=locVLA))
                    self.az_i = altaz.az.deg
                    self.alt_i = altaz.alt.deg
                    hrang = (15*self.lst) - self.obs_ra  # deg
                    if hrang > 180:
                        hrang -= 360
                    if hrang < -180:
                        hrang += 360
                    tmp1 = np.sin(hrang*DEG2RAD)
                    tmp2 = np.tan(VLALAT*DEG2RAD)*np.cos(self.obs_dec*DEG2RAD) - \
                        np.sin(self.obs_dec*DEG2RAD)*np.cos(hrang*DEG2RAD)
                    self.parang_i = np.arctan2(tmp1, tmp2)*RAD2DEG
                    if self.duration is not None:
                        t_end = Time(self.mjdtime+(self.duration /
                                                   86400.), format='mjd')  # end time
                        lst_end = t_end.sidereal_time(
                            'apparent', VLALON).hour  # hrs
                        altaz = coord.transform_to(
                            AltAz(obstime=t_end, location=locVLA))
                        self.az_f = altaz.az.deg
                        self.alt_f = altaz.alt.deg
                        hrang = (15*lst_end) - self.obs_ra  # deg
                        if hrang > 180:
                            hrang -= 360
                        if hrang < -180:
                            hrang += 360
                        tmp1 = np.sin(hrang*DEG2RAD)
                        tmp2 = np.tan(VLALAT*DEG2RAD)*np.cos(self.obs_dec*DEG2RAD) - \
                            np.sin(self.obs_dec*DEG2RAD)*np.cos(hrang*DEG2RAD)
                        self.parang_f = np.arctan2(tmp1, tmp2)*RAD2DEG
        except KeyError:
            self.mjdtime = None
            self.lst = None
            self.az_i = None
            self.alt_i = None
            self.parang_i = None
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
        """Sets image cycle, semester, and ass_flag"""
        self.ass_flag = False
        config = self.config
        if config is None:
            self.semester = None
            self.cycle = None
        else:
            for cycle in res_dict[config].keys():
                if self.mjdtime <= res_dict[config][cycle]['mjd'][1] and self.mjdtime > res_dict[config][cycle]['mjd'][0]:
                    self.semester = res_dict[config][cycle]['semester']
                    self.cycle = cycle
                    if self.bmin <= res_dict[config][cycle]['bmin'][1] and self.bmin >= res_dict[config][cycle]['bmin'][0]:
                        self.ass_flag = True
                    break
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
    id : int
        Uniquely identifies the source in the **assoc_source** table.
    res_class : str
        Resolution class of the image in which the source was found.
        'A': res <= 15", 'B': 15" < res <= 35", 'C': 35" < res <= 60",
        'D': res > 60".
    ndetect : int
        Number of times this same source has been detected in other
        VLITE images.
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
        # assoc_source attributes
        self.id = None
        self.res_class = None
        self.ndetect = None
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
        from the image pointing center in degrees and sets
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
        world = np.array([[self.ra, self.dec]])
        pcoordS = imobj.wcsobj.wcs_world2pix([[self.ra, self.dec]], 1)
        world = np.array([[imobj.obs_ra, imobj.obs_dec]])
        pcoordI = imobj.wcsobj.wcs_world2pix(world, 1)
        theta = np.arctan2(pcoordI[0][0]-pcoordS[0]
                           [0], pcoordS[0][1]-pcoordI[0][1])  # rad
        self.polar_angle = theta * 180.0/np.pi  # deg

    def correct_flux(self, pri_freq):
        """Applies a 1-D radial correction to all flux measurements
        from PyBDSF. The correction factor has been empirically derived
        from flux comparisons of sources as a function of angular
        distance from the image center. The scale factor increases
        farther from the pointing center as the beam power decreases.

        Parameters
        ----------
        pri_freq : float
            Primary frequency of the observations in GHz.
        """
        # Correct for beam response - only 1D (sym. beam) for now
        pb_power, pb_err = beam_tools.find_nearest_pbcorr(self.dist_from_center,
                                                          pri_freq)
        # List of attributes to correct
        attrs = ['total_flux', 'peak_flux', 'total_flux_isl', 'rms_isl', 'mean_isl',
                 'resid_rms', 'resid_mean']
        for attr in attrs:
            corr_val = getattr(self, attr) / pb_power
            setattr(self, attr, corr_val)
        # Add systematic uncertainty to all flux errors
        self.e_total_flux = np.sqrt(
            (self.total_flux * pb_err)**2. + self.e_total_flux**2.)
        self.e_peak_flux = np.sqrt(
            (self.peak_flux * pb_err)**2. + self.e_peak_flux**2.)
        self.total_flux_islE = np.sqrt(
            (self.total_flux_isl * pb_err)**2. + self.total_flux_islE**2.)

    def calc_snr(self):
        # Compute S/N of source detection
        self.snr = (self.peak_flux - self.mean_isl) / self.rms_isl

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
        attrs = ['total_flux', 'peak_flux']
        for attr in attrs:
            corr_val = getattr(self, attr) / pb_power
            setattr(self, attr, corr_val)
        # Compute S/N of source detection
        #self.snr = (self.peak_flux - self.mean_isl) / self.rms_isl
        # Add systematic uncertainty to all flux errors
        self.e_total_flux = np.sqrt(
            (self.total_flux * pb_err)**2. + self.e_total_flux**2.)
        self.e_peak_flux = np.sqrt(
            (self.peak_flux * pb_err)**2. + self.e_peak_flux**2.)
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
        if imobj.config == 'A' or imobj.config == 'BnA':
            a0 = 1.054294
            a1 = 4.497675e-3
            c0 = 18.277850
            c1 = -1.917682
        elif imobj.config == 'B' or imobj.config == 'CnB':
            a0 = 1.028851
            a1 = 2.17631e-3
            c0 = 10.493944
            c1 = -1.804784
        elif imobj.config == 'C' or imobj.config == 'DnC':
            a0 = 1.007876
            a1 = 1.761468e-3
            c0 = 12.600512
            c1 = -2.002195
        elif imobj.config == 'D':  # need to determine D config values
            a0 = 1.0
            a1 = 0.
            c0 = 0.
            c1 = 0.
        else:  # unknown config
            a0 = 1.0
            a1 = 0.
            c0 = 0.
            c1 = 0.
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
    for oldsrc in out.sources:
        newsrcs.append(DetectedSource())
        newsrcs[-1].cast(oldsrc)
        newsrcs[-1].image_id = img.id
        newsrcs[-1].calc_center_dist(img)
        newsrcs[-1].calc_polar_angle(img)

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


def set_nsn(img):
    """Counts number of SN tables in UVOUT file

    Parameters
    ----------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute
        values set from header info.

    Returns
    -------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute
        nsn updated.
    """
    # Skip for VCSS images
    if 'VCSS' in img.filename or 'vcss' in img.filename:
        img.nsn = None
        return img

    # Determine dir with uvout file
    a = img.filename.split('/')
    uvdir = ''
    for i in range(len(a)-2):
        uvdir += (a[i]+'/')
    uvdir += 'UVFiles/'
    # Check dir
    if not os.path.isdir(uvdir):
        img.nsn = None
        dbclasses_logger.info('Did not find UVFiles dir {}.'.format(uvdir))
        return img

    # Determine name of uvout file
    b = a[-1].split('.')
    uvname1 = ''
    for i in range(len(b)-2):
        uvname1 += (b[i]+'.')
    uvname = uvname1 + 'uvout.gz'
    # Check file
    if not os.path.isfile(os.path.join(uvdir, uvname)):
        # unzipped?
        uvname = uvname1 + 'uvout'
        if not os.path.isfile(os.path.join(uvdir, uvname)):
            img.nsn = None
            dbclasses_logger.info(
                'Did not find uvout file, with or w/o .gz {}'.format(uvname))
            return img

    # Count number of SN tables
    cnt = 0
    with fits.open(os.path.join(uvdir, uvname)) as hdu:
        for i in range(1, len(hdu)):
            if hdu[i].header['EXTNAME'] == 'AIPS SN':
                cnt += 1
    img.nsn = cnt

    dbclasses_logger.info('{}: nsn= {}'.format(uvname, img.nsn))

    return img


def init_image(impath, alwaysass):
    """Initializes an object of the Image class and sets values 
    for its attributes from the fits file header using
    the ``header_attrs`` object method.

    Parameters
    ----------
    alwaysass : 'always associate' option. If True image sources
              will always be associated
    """
    imobj = Image()
    imobj.process_image(impath)
    hdu, hdr = imobj.read()
    # Use header info to set attributes
    imobj.header_attrs(hdr)
    # Fix header keywords for PyBDSF, if necessary
    header_changed = False
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

    # Set nsn
    #dbclasses_logger.info('calling set_nsn {}'.format(imobj.filename))
    imobj = set_nsn(imobj)

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
