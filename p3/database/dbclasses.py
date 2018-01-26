"""This module defines classes and their methods for
Image and DetectedSource objects.

"""
import os
import re
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import ephem
import pybdsfcat
from sourcefinding import beam_tools


class Image(object):
    """Creates an Image object that can be easily 
    passed to SQL table insertion functions. Object 
    attributes correspond to the database image table 
    column values.

    Parameters
    ----------
    image : str
        Directory path to the fits image file.

    Attributes
    ----------
    id : int
        Numerical index. Incremented by `PostgreSQL` upon insertion.
    filename : str
        Full path name for the image fits file.
    imsize : str
        Image size in pixels -- (NAXIS1, NAXIS2).
    obs_ra : float
        Right ascension of image pointing center in degrees.
    obs_dec : float
        Declination of image pointing center in degrees.
    pixel_scale : float
        Spatial conversion from pixel to angular size.
        Units are arcseconds / pixel.
    obj : str
        Name of observed object.
    obs_date : str
        Date when observations were taken.
    map_date : str
        Date when data was imaged.
    obs_freq : float
        Frequency of the observations in MHz (actual
        frequency of the image).
    pri_freq : float
        Frequency of the simultaneously acquired primary
        band data in GHz (NOT the frequency of the image
        on hand).
    bmaj : float
        Size of the beam semi-major axis in arcseconds.
    bmin : float
        Size of the beam semi-minor axis in arcseconds.
    bpa : float
        Beam position angle in degrees.
    noise : float
        Estimate of the image noise in mJy/beam.
    peak : float
        Peak brightness in the image in mJy/beam.
    config : str
        Very Large Array configuration.
    nvis : int
        Number of visibilities in the data before imaging.
    mjdtime : float
        Modified Julian Date (days since 0h Nov 17, 1858).
    tau_time : float
        Integration time on source in seconds.
    duration : float
        Total time spanned by the observations in seconds.
    nsrc : int
        Number of sources extracted during source finding.
    rms_box : str
        `PyBDSF` parameter used during source finding.
        From `PyBDSF` docs: "The first integer, boxsize, is the 
        size of the 2-D sliding box [in pixels] for calculating 
        the rms and mean over the entire image. The second, stepsize, 
        is the number of pixels by which this box is moved for the next 
        measurement."
    error_id : int
        Numerical code assigned when an image fails a quality
        check. Each number has a corresponding explanation in
        the database error table.
    stage : int
        Highest processing stage completed. 1 = Image object initialized,
        2 = source finding, 3 = source association, 4 = catalog matching.
    radius : float
        Radius defining the circular field-of-view in degrees in which
        sources were extracted.
    """
    # A class variable to count the number of images
    num_images = 0
    
    def __init__(self, image):
        pri = re.findall('\/([0-9.]+[A-Z])', image)
        if pri[0][-1:] == 'M':
            pri_freq = float(pri[0][:-1]) / 1000. # GHz
        else:
            pri_freq = float(pri[0][:-1]) # GHz
        self.id = None
        self.filename = image
        self.imsize = None
        self.obs_ra = None
        self.obs_dec = None
        self.pixel_scale = None
        self.obj = None
        self.obs_date = None
        self.map_date = None
        self.obs_freq = None
        self.pri_freq = pri_freq
        self.bmaj = None
        self.bmin = None
        self.bpa = None
        self.noise = None
        self.peak = None
        self.config = None
        self.nvis = None
        self.mjdtime = None
        self.tau_time = None
        self.duration = None
        self.nsrc = None
        self.rms_box = None
        self.error_id = None
        self.stage = 1
        self.radius = None
        self.nearest_problem = None
        self.separation = None

        # Increase the image count by one
        Image.num_images += 1


    @classmethod
    def image_count(cls):
        """Returns the number of images initialized."""
        return cls.num_images


    def read(self):
        """Reads fits image data and header."""
        try:
            data, header = fits.getdata(self.filename, header=True)
            return data, header
        except:
            print('\nERROR: Problem reading image.\n')


    def write(self, data, header, owrite=False):
        """Writes new fits image data and header."""
        fits.writeto(self.filename, data, header, overwrite=owrite)


    def image_qa(self, params):
        """Performs quality checks on the image pre-source finding
        to flag images with the following issues:
        1.) integration time on source (tau_time) < 60 s
        2.) noise < 0 mJy/beam or > 1000 mJy/beam
        3.) beam semi-major/semi-minor axis > 3 (too elliptical)
        4.) pointing center within 20 degrees of a known bright
        radio source (Cas A, Cygnus A, Taurus A, Hercules A,
        Virgo A (M87), the Sun, Jupiter, Galactic Center)

        """
        # QA requirements
        min_tos = params['min time on source (s)']
        max_rms = params['max noise (mJy/beam)']
        max_ellip = params['max beam axis ratio']
        prob_sep = params['min problem source separation (deg)']
        
        # VLITE can't track planets, so flag all of them
        sun = ephem.Sun()
        mer = ephem.Mercury()
        ven = ephem.Venus()
        mars = ephem.Mars()
        jup = ephem.Jupiter()
        sat = ephem.Saturn()
        ur = ephem.Uranus()
        nep = ephem.Neptune()
        plu = ephem.Pluto()
        # convert from MJD to JD
        jd = self.mjdtime + 2400000.5
        jdt = Time(jd, format='jd', scale='utc')
        sun.compute(jdt.iso)
        mer.compute(jdt.iso)
        ven.compute(jdt.iso)
        mars.compute(jdt.iso)
        jup.compute(jdt.iso)
        sat.compute(jdt.iso)
        ur.compute(jdt.iso)
        nep.compute(jdt.iso)
        plu.compute(jdt.iso)
        bright_sources = {'Sun' : SkyCoord(str(sun.a_ra), str(sun.a_dec),
                                           unit=('hourangle', 'deg')),
                          'Mercury' : SkyCoord(str(mer.a_ra), str(mer.a_dec),
                                               unit=('hourangle', 'deg')),
                          'Venus' : SkyCoord(str(ven.a_ra), str(ven.a_dec),
                                             unit=('hourangle', 'deg')),
                          'Mars' : SkyCoord(str(mars.a_ra), str(mars.a_dec),
                                            unit=('hourangle', 'deg')),
                          'Jupiter' : SkyCoord(str(jup.a_ra), str(jup.a_dec),
                                               unit=('hourangle', 'deg')),
                          'Saturn' : SkyCoord(str(sat.a_ra), str(sat.a_dec),
                                              unit=('hourangle', 'deg')),
                          'Uranus' : SkyCoord(str(ur.a_ra), str(ur.a_dec),
                                              unit=('hourangle', 'deg')),
                          'Neptune' : SkyCoord(str(nep.a_ra), str(nep.a_dec),
                                               unit=('hourangle', 'deg')),
                          'Pluto' : SkyCoord(str(plu.a_ra), str(plu.a_dec),
                                             unit=('hourangle', 'deg')),
                          'Cas A' : SkyCoord(350.866250, 58.811667, unit='deg'),
                          'Cyg A' : SkyCoord(299.867917, 40.733889, unit='deg'),
                          'Tau A' : SkyCoord(83.633333, 22.014444, unit='deg'),
                          'Her A' : SkyCoord(252.783750, 4.992500, unit='deg'),
                          'M87' : SkyCoord(187.705833, 12.391111, unit='deg'),
                          'GC' : SkyCoord(266.416833, -29.007806, unit='deg')}

        image_center = SkyCoord(self.obs_ra, self.obs_dec, unit='deg')

        print('\nPerforming preliminary image quality checks...')

        if self.tau_time is not None:
            if self.tau_time < min_tos:
                print('IMAGE FAILED QA: integration time on source < {} s'.
                      format(min_tos))
                self.error_id = 1
                return
            else: pass
        else:
            if self.duration is not None and self.duration < min_tos:
                print('IMAGE FAILED QA: integration time on source < {} s'.
                      format(min_tos))
                self.error_id = 1
                return
            else: pass

        if self.noise is not None:
            if self.noise < 0. or self.noise > max_rms:
                print('IMAGE FAILED QA: image noise is {}. Max allowed is {}'.
                      format(self.noise, max_rms)) 
                self.error_id = 2
                return
            else: pass
        else:
            pass

        try:
            axis_ratio = self.bmaj / self.bmin
        except TypeError:
            axis_ratio = 1.
        if axis_ratio > max_ellip:
            print('IMAGE FAILED QA: beam axis ratio is {}. Max allowed is {}'.
                  format(axis_ratio, max_ellip))
            self.error_id = 3
            return
        else: pass

        min_sep = 999999.99
        for src, loc in bright_sources.items():
            ang_sep = image_center.separation(loc).degree
            if ang_sep < min_sep:
                min_sep = ang_sep
                self.separation = ang_sep
                self.nearest_problem = src
            else: pass
            if ang_sep < prob_sep:
                print('IMAGE FAILED QA: {} is within {} degrees'.format(
                    src, prob_sep))
                self.error_id = 4
                return
            else: pass
        # Check to see if image header says it's looking at the NCP
        if self.obj == 'NCP':
            print('IMAGE FAILED QA: pointing is at NCP')
            self.error_id = 4
            return
        else: pass
        # Some planet observations have obs_ra, obs_dec = 0, 0
        if self.obs_ra == 0. and self.obs_dec == 0.:
            print('IMAGE FAILED QA: pointing center is at RA, Dec = 0, 0')
            self.error_id = 4
            return
        else: pass

        print('...image passed.')


    def source_qa(self, sources, params):
        """Approximates the expected number of sources in an
        image based on source counts from WENSS scaled by
        noise and fitted beam. See EP's PythonTools.

        """
        max_src_metric = params['max source metric']
        
        print('\nPerforming source count quality checks...')
        # First, check if there are any sources
        if len(sources) < 1:
            print('IMAGE FAILED QA: 0 sources were detected')
            self.error_id = 6
            return
        else: pass

        # Get number of actual sources within central 1.5 degrees
        image_center = SkyCoord(self.obs_ra, self.obs_dec, unit='deg')
        nsrc_cut = 0
        for src in sources:
            src_loc = SkyCoord(src.ra, src.dec, unit='deg')
            deg_from_center = image_center.separation(src_loc).degree
            if deg_from_center <= 1.5:
                nsrc_cut += 1

        # Estimate number of expected sources within central 1.5 degrees
        nsrc_exp = beam_tools.expected_nsrc(self.noise)
        # Compute metric for determining if source count is way off
        nsrc_metric = (float(nsrc_cut) - nsrc_exp) / nsrc_exp
        
        if nsrc_metric > max_src_metric:
            print('IMAGE FAILED QA: too many sources detected, ')
            print('nsrc_metric = {}; limit = {}'.format(
                nsrc_metric, max_src_metric))
            self.error_id = 7
            return
        else: pass

        print('...image passed.')

        
    def header_attrs(self, hdr):
        """This method does all the heavy lifting of 
        extracting and assigning the object attributes 
        from the fits image header and assigns ``None`` 
        if the keyword is missing in the header.

        """
        try:
            naxis1 = hdr['NAXIS1']
            naxis2 = hdr['NAXIS2']
            self.imsize = str((naxis1, naxis2)) # pixels
        except KeyError:
            self.imsize = None
        try:
            self.obs_ra = hdr['OBSRA'] # deg
        except KeyError:
            self.obs_ra = None
        try:
            self.obs_dec = hdr['OBSDEC'] # deg
        except KeyError:
            self.obs_dec = None
        try:
            self.pixel_scale = abs(hdr['CDELT1']) * 3600. # arcsec/pixel
        except KeyError:
            try:
                self.pixel_scale = abs(hdr['CDELT2']) * 3600. # arcsec/pixel
            except KeyError:
                self.pixel_scale = None
        try:
            self.obj = hdr['OBJECT']
        except KeyError:
            self.obj = None
        try:
            self.obs_date = hdr['DATE-OBS']
        except KeyError:
            self.obs_date = None
        try:
            self.map_date = hdr['DATE-MAP']
        except KeyError:
            self.map_date = None
        try:
            self.obs_freq = hdr['RESTFREQ'] / 10**6. # MHz
        except KeyError:
            try:
                if hdr['CTYPE3'] == 'FREQ':
                    self.obs_freq = hdr['CRVAL3'] / 10**6. # MHz
                else:
                    self.obs_freq = hdr['CRVAL4'] / 10**6. # MHz
            except KeyError:
                self.obs_freq = None
        try:
            self.bmaj = hdr['BMAJ'] * 3600. # arcsec
            self.bmin = hdr['BMIN'] * 3600. # arcsec
            self.bpa = hdr['BPA'] # deg
        except KeyError:
            try:
                self.bmaj = hdr['CLEANBMJ'] * 3600. # arcsec
                self.bmin = hdr['CLEANBMN'] * 3600. # arcsec
                self.bpa = hdr['CLEANBPA'] # deg
            except KeyError:
                try:
                    # Search for beam params in AIPS history
                    hl = list(hdr['HISTORY'])
                    for line in hl:
                        x = re.findall('BMAJ=\s+([0-9]\S+)', line)
                        y = re.findall('BMIN=\s+([0-9]\S+)', line)
                        z = re.findall('BPA=\s+([0-9]\S+)', line)
                        if len(x) > 0:
                            self.bmaj = float(x[0]) * 3600. # arcsec
                        if len(y) > 0:
                            self.bmin = float(y[0]) * 3600. # arcsec
                        if len(z) > 0:
                            self.bpa = float(z[0]) # deg
                except KeyError:
                    self.bmaj = None
                    self.bmin = None
                    self.bpa = None
        try:
            self.noise = hdr['ACTNOISE'] * 1000. # mJy/beam
        except KeyError:
            self.noise = None
        try:
            self.peak = hdr['PEAK'] * 1000. # mJy/beam
        except KeyError:
            try:
                self.peak = hdr['DATAMAX'] * 1000. # mJy/beam
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
        try:
            self.mjdtime = hdr['MJDTIME']
        except KeyError:
            self.mjdtime = None
        try:
            self.tau_time = hdr['TAU_TIME'] # sec
        except KeyError:
            self.tau_time = None
        try:
            self.duration = hdr['DURATION'] # sec
        except KeyError:
            self.duration = None


    def log_attrs(self, pybdsfdir):
        """Image class method to assign values to object 
        attributes which come from source finding using info 
        in the `PyBDSF` output log file. This is useful when 
        source finding has already been run on an image and the 
        results are being recorded after the fact using the 
        `PyBDSF` output files.

        Parameters
        ----------
        pybdsfdir : str
            Directory path to the `PyBDSF` output files.

        Returns
        -------
        rms_box : str
            `PyBDSF` parameter used during source finding.
            From `PyBDSF` docs: "The first integer, boxsize, is the 
            size of the 2-D sliding box [in pixels] for calculating 
            the rms and mean over the entire image. The second, stepsize, 
            is the number of pixels by which this box is moved for the next 
            measurement."
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
        gresid_std = float(y[-1]) * 1000. # mJy/beam
        z = re.findall('.*raw rms =\s+(.*) mJy', log)
        raw_rms = float(z[-1])
        return rms_box, gresid_std, raw_rms


class DetectedSource(object):
    """Class of objects to store properties of sources extracted
    from an image by `PyBDSF`. Attributes are translated from
    attributes of the class of objects output by `PyBDSF`.

    References
    ----------
    Refer to the `PyBDSF` documentation ([1]_) for definitions 
    of their output columns.

    .. [1] http://www.astron.nl/citt/pybdsm/write_catalog.html#definition-of-output-columns
    """
    def __init__(self):
        self.src_id = None
        self.isl_id= None
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
        # assoc_source attributes
        self.id = None
        self.beam = None
        self.ndetect = None
        self.nmatches = None
        # vlite_unique attribute
        self.detected = None


    def cast(self, origsrc):
        """Attributes of a `bdsf.gaul2srl.Source` object
        output from `bdsf.process_image` are cast to define
        attributes of a DetectedSource object.

        Parameters
        ----------
        origsrc : bdsf.gaul2srl.Source instance
            Object which stores all the measured properties
            of a source extracted from an image in `PyBDSF`.

        """
        self.src_id = origsrc.source_id
        self.isl_id= origsrc.island_id
        self.ra = origsrc.posn_sky_centroid[0] # deg
        self.e_ra = origsrc.posn_sky_centroidE[0] # deg
        self.dec = origsrc.posn_sky_centroid[1] # deg
        self.e_dec = origsrc.posn_sky_centroidE[1] # deg
        self.total_flux = origsrc.total_flux * 1000. # mJy
        self.e_total_flux = origsrc.total_fluxE * 1000. # mJy
        self.peak_flux = origsrc.peak_flux_centroid * 1000. # mJy/beam
        self.e_peak_flux = origsrc.peak_flux_centroidE * 1000. # mJy/beam
        self.ra_max = origsrc.posn_sky_max[0] # deg
        self.e_ra_max = origsrc.posn_sky_maxE[0] # deg
        self.dec_max = origsrc.posn_sky_max[1] # deg
        self.e_dec_max = origsrc.posn_sky_maxE[1] # deg
        self.maj = origsrc.size_sky[0] * 3600. # arcsec
        self.e_maj = origsrc.size_skyE[0] * 3600. # arcsec
        self.min = origsrc.size_sky[1] * 3600. # arcsec
        self.e_min = origsrc.size_skyE[1] * 3600. # arcsec
        self.pa = origsrc.size_sky[2] # deg
        self.e_pa = origsrc.size_skyE[2] # deg
        self.dc_maj = origsrc.deconv_size_sky[0] * 3600. # deconvolved
        self.e_dc_maj = origsrc.deconv_size_skyE[0] * 3600.
        self.dc_min = origsrc.deconv_size_sky[1] * 3600.
        self.e_dc_min = origsrc.deconv_size_skyE[1] * 3600.
        self.dc_pa = origsrc.deconv_size_sky[2]
        self.e_dc_pa = origsrc.deconv_size_skyE[2]
        self.total_flux_isl = origsrc.total_flux_isl * 1000. # mJy
        self.total_flux_islE = origsrc.total_flux_islE * 1000. # mJy
        self.rms_isl = origsrc.rms_isl * 1000. # mJy/beam
        self.mean_isl = origsrc.mean_isl * 1000. # mJy/beam
        self.resid_rms = origsrc.gresid_rms * 1000. # mJy/beam
        self.resid_mean = origsrc.gresid_mean * 1000. # mJy/beam
        self.code = origsrc.code


    def correct_flux(self, imobj):
        img_center = SkyCoord(imobj.obs_ra, imobj.obs_dec, unit='deg')
        src_loc = SkyCoord(self.ra, self.dec, unit='deg')
        angle = img_center.separation(src_loc).degree
        # Correct for beam response - only 1D (sym. beam) for now
        pb_power = beam_tools.find_nearest_pbcorr(angle)
        # List of attributes to correct
        attrs = ['total_flux', 'e_total_flux', 'peak_flux', 'e_peak_flux',
                 'total_flux_isl', 'total_flux_islE', 'rms_isl', 'mean_isl',
                 'resid_rms', 'resid_mean']
        for attr in attrs:
            corr_val = getattr(self, attr) / pb_power
            setattr(self, attr, corr_val)


def dict2attr(obj, dictionary):
    """Sets dictionary key, value pairs to object attributes.
    
    Parameters
    ----------
    obj : class instance
        Any object initialized from a class.
    dictionary : dict
        Dictionary key, value pairs become attribute, 
        value pairs.
    """
    for key in dictionary.keys():
        setattr(obj, key, dictionary[key])


def translate_from_txtfile(img, pybdsfdir):
    """Creates a list of DetectedSource objects
    from reading in a source list text file output
    by `PyBDSF`.
    
    Parameters
    ----------
    img : database.dbclasses.Image instance
        Initialized Image object with attribute
        values set from header info.
    pybdsfdir : str
        Directory path to `PyBDSF` files.

    Returns
    -------
    img : database.dbclasses.Image instance
        Initialized Image object with attribute
        values updated with source finding results.
    sources : list of DetectedSource objects
        DetectedSource objects with attribute values
        set from `PyBDSF` output source list text file.
    """
    # Set the file path
    prefix = re.findall('.*\/(.*)\.', img.filename)[0]
    srl = prefix + '.pybdsm.srl'
    catalog = os.path.join(pybdsfdir, srl)
    # Read the catalog
    print('\nExtracting sources from {}'.format(srl))
    sources = pybdsfcat.read_catalog(catalog)

    # Count the number of sources
    img.nsrc = len(sources)

    # Extract rms_box size and noise from log file
    img.rms_box, img.noise, raw_rms = img.log_attrs(pybdsfdir)

    return img, sources


def translate(img, out):
    """Method to translate `PyBDSF` output within the
    pipeline to DetectedSource objects and update 
    the Image object.

    Parameters
    ----------
    img : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info.
    out : bdsf.image.Image instance
        The object output by `PyBDSF` after running its
        source finding task `process_image()`. Contains
        a list of bdsf.gaul2srl.Source objects which
        are translated into DetectedSource objects.

    Returns
    -------
    img : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    newsrcs : list
        List of database.dbclasses.DetectedSource objects.
        Attributes of each object are set from the `PyBDSF`
        output object.    
    """
    # Add PyBDSF defined attributes
    img.nsrc = out.nsrc
    img.rms_box = str(out.rms_box)
    # Try updating any missing attributes from header info
    # using PyBDSF's output object
    if img.imsize is None:
        img.imsize = str(out._original_shape) # pixels
    if img.obs_freq is None:
        img.obs_freq = out.frequency / 10**6. # MHz
    if img.bmaj is None:
        img.bmaj = out.beam[0] * 3600. # arcsec
    if img.bmin is None:
        img.bmin = out.beam[1] * 3600. # arcsec
    if img.bpa is None:
        img.bpa = out.beam[2] # deg
    if img.noise is None:
        img.noise = np.std(out.resid_gaus_arr) * 1000. # mJy/beam

    # Translate PyBDSF output source objects to our own
    # DetectedSource objects
    newsrcs = []
    for oldsrc in out.sources:
        newsrcs.append(DetectedSource())
        newsrcs[-1].cast(oldsrc)
        newsrcs[-1].image_id = img.id

    return img, newsrcs
