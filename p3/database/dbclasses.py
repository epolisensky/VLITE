import sys
import os
import re
from astropy.io import fits
from astropy.coordinates import Angle
import astropy.units as u


class Image(object):
    """Creates a database Image object that can be easily 
    passed to SQL table insertion functions. Object attributes
    correspond to the Image table column values."""

    # A class variable to count the number of images
    num_images = 0
    
    def __init__(self, image):
        self.filename = image
        self.imsize = None
        self.obs_ra = None
        self.obs_dec = None
        self.pixel_scale = None
        self.obj = None
        self.obs_date = None
        self.map_date = None
        self.freq = None
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

        # Increase the image count by one
        Image.num_images += 1


    @classmethod
    def image_count(cls):
        """Prints the number of images initialized."""
        print('\nProcessed {:d} images.'.format(cls.num_images))


    def read(self):
        """Reads fits image data and header."""
        try:
            data, header = fits.getdata(self.filename, header=True)
            return data, header
        except:
            print('\nERROR: Problem reading image.\n')


    def write(self, data, header):
        """Writes new fits image data and header."""
        try:
            fits.writeto(self.filename, data, header, overwrite=True)           
        except:
            print('\nNo image changes to write.')


    def search_radius(self):
        try:
            imsize = float(self.imsize[1:].split(',')[0])
            pixscale = self.pixel_scale / 3600. # degrees
            radius = Angle(imsize * pixscale * u.deg) / 2.
            return radius.degree
        except: # Could be TypeError or AttributeError
            data, hdr = self.read()
            radius = Angle(hdr['NAXIS1'] * hdr['CDELT2'] * u.deg) / 2.
            return radius.degree


    def header_attrs(self, hdr):
        """This method does all the heavy lifting of extracting and
        assigning the desired object attributes (table column values)
        from the fits image header and assigns 'None' if the keyword
        is missing in the header."""
        
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
            self.freq = hdr['RESTFREQ'] / 10**6. # MHz
        except KeyError:
            try:
                if hdr['CTYPE3'] == 'FREQ':
                    self.freq = hdr['CRVAL3'] / 10**6. # MHz
                else:
                    self.freq = hdr['CRVAL4'] / 10**6. # MHz
            except KeyError:
                self.freq = None
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
    """This is basically just a convenience class created
    to allow consistent naming and access of properties of
    sources detected by PyBDSF across all data input types 
    (i.e. inside from a bdsf.process_image output object or 
    outside from a database table or text file)."""
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
        self.catalog_id = None
        self.match_id = None
        self.min_deRuiter = None
        #self.sig = None
	#self.alpha = None
	#self.phi = None


    def cast(self, origsrc):
        """Attributes of a bdsf.gaul2srl.Source object
        output from bdsf.process_image are 'cast' to define
        attributes of a DetectedSource object."""
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
        #self.sig = self.peak_flux / self.rms_isl
