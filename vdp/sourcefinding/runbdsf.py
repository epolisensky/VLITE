"""This module contains all the functionality to run a FITS image
through the LOFAR source finding software PyBDSF 
(https://github.com/lofar-astron/PyBDSF). The class
``BDSFImage()`` creates an object whose attributes can be any number
of PyBDSF parameters. Source finding is performed by calling the
object method ``find_sources()``, which calls the PyBDSF function
``process_image()`` and returns the output object.

"""
import warnings
import logging
from datetime import datetime
import numpy as np
import bdsf
from database.dbclasses import Image
from timeout import timeout


# create logger
sf_logger = logging.getLogger('vdp.sourcefinding.runbdsf')


def write_sources(out):
    """Calls PyBDSF functions ``write_catalog()``
    and ``export_image()`` to record results from 
    source finding on an image.
    
    Parameters
    ----------
    out : ``bdsf.image.Image`` instance
        The object output by PyBDSF after running its
        source finding task ``process_image()``.
    """
    # Write the source list catalog, ascii and ds9 regions file
    out.write_catalog(format='ds9', catalog_type='srl', clobber=True)
    #out.write_catalog(format='ascii', catalog_type='srl', clobber=True)

    # Write the residual image
    #out.export_image(img_type='gaus_resid', clobber=True)
    # Write the model image
    #out.export_image(img_type='gaus_model', clobber=True)


class BDSFImage(Image):
    """Object to be manipulated and read into PyBDSF.
    Inherits all methods defined for Image class, but
    overrides initialization. The number of attributes
    and their values will change based on what the user
    specifies in the run configuration file.

    """
    def __init__(self, image, **kwargs):
        """Initializes the Image subclass object. PyBDSF will 
        use any arguments it recognizes and ignore the rest."""
        self.filename = image
        # These will be overwritten if in config file
        self.quiet = True
        self.box_incr = 20
        self.max_iter = 5
        # Set our own default rms_box parameter
        self.set_rms_box()
        # Setting force-fitting coords and stop_at to None 
        #  will not affect blind source fitting
        self.src_ra_dec = None
        self.stop_at = None
        # Set attributes from config file
        for key, value in kwargs.items():
            if key == 'rms_box':
                if value == '' or value == 'None':
                    value = None
            if key == 'rms_box_bright':
                if value == '' or value == 'None':
                    value = None
            if key == 'adaptive_thresh':
                if value == '' or value == 'None':
                    value = None
            setattr(self, key, value)


    def set_rms_box(self):
        """Sets the PyBDSF ``rms_box`` parameter to a box size
        1/10th of the image size and a step size one third of
        the box size. This "VLITE default" ``rms_box`` yields slightly
        better results with fewer artifacts than if left as a
        free parameter for PyBDSF to calculate. A custom
        ``rms_box`` can be defined in the configuratrion file
        which will supersede the one defined here.

        """
        data, hdr = Image.read(self)
        box_size = int(round(hdr['NAXIS2'] / 10.))
        step_size = int(round(box_size / 3.))
        self.rms_box = (box_size, step_size)
        small_box_size = int(round(box_size / 5.))
        small_step_size = int(round(step_size / 5.))
        self.rms_box_bright = (small_box_size, small_step_size)
      

    def get_attr(self):
        """Return all object attributes as a dictionary.
        This is fed into ``bdsf.process_image()`` as parameter
        arguments.

        """
        return self.__dict__

    # Function will timeout after 5 min
    @timeout()
    def find_sources(self):
        """Run PyBDSF ``process_image()`` task using object 
        attributes as parameter inputs. Returns ``None`` if 
        PyBDSF fails. Wrapped in a timeout function so
        processing is killed if taking longer than 5 minutes.

        """
        sf_logger.info('Extracting sources...')
        start = datetime.now()
        opts = self.get_attr()
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', r'invalid value')
            try:
                out = bdsf.process_image(opts)
            except:
                out = None

        try:
            sf_logger.info(' -- found {} sources in {:.2f} seconds'.format(
                out.nsrc, (datetime.now() - start).total_seconds()))
        except AttributeError:
            try:
                sf_logger.info(' -- found {} islands in {:.2f} seconds'.format(
                    out.nisl, (datetime.now() - start).total_seconds()))
            except AttributeError:
                sf_logger.info(' -- PyBDSF failed to process image.')

        return out
            

    def minimize_islands(self):
        """Incrementally increases and decreases the ``bdsf.process_image()``
        argument ``rms_box`` box size until a minimum number of islands are 
        found. This technique works best on images with significant 
        artifacts where false detections are the biggest concern.
        This also takes a really long time since it is running
        ``bdsf.process_image`` numerous times, so is probably only
        ever worth using when analyzing a small number of images.

        """
        # Initial run
        sf_logger.info('Starting minimize_islands with box {}...'.format(
            self.rms_box))
        self.stop_at = 'isl' # stop fitting at islands
        out = self.find_sources()
        if out is not None:
            box0 = out.rms_box # record initial box
            min_isl = out.nisl # initialize island number minimum
        else:
            # Only fails when user's box is too small
            box0 = self.rms_box            
            min_isl = 99999
            sf_logger.info('PyBDSF has failed with box {}.'.format(box0))
            sf_logger.info('Increasing rms_box...')

        opt_box = box0 # initialize optimal box
        box_size = box0[0]

        # Begin increasing box size loop
        i = 0
        while i < self.max_iter: # stop at max number of iterations
            box_size += self.box_incr
            box_step = int(box_size / 3.)
            self.rms_box = (box_size, box_step)
            sf_logger.info('Trying box {}...'.format(self.rms_box))
            out = self.find_sources()
            if out is not None:
                if out.nisl < min_isl: # new winner, keep going
                    min_isl = out.nisl
                    opt_box = self.rms_box
                elif out.nisl == min_isl: # tie, keep going
                    pass
                else: # number of islands increased, stop
                    break
            else:
                # usually happens when box size is too small and "an
                # unphysical rms value was encountered", so keep going
                sf_logger.info('PyBDSF has failed with box {}.'.format(
                    self.rms_box))
                sf_logger.info('Increasing rms_box...')
            i += 1

        # Begin decreasing box size loop
        i = 0
        box_size = box0[0] # reset the box size back to original
        while i < self.max_iter:
            box_size -= self.box_incr
            box_step = int(box_size / 3.)
            self.rms_box = (box_size, box_step)
            if box_size > 0: # stop if box size goes to 0 or below
                sf_logger.info('Trying box {}...'.format(self.rms_box))
                out = self.find_sources()
            else:
                break
            if out is not None:
                if out.nisl < min_isl: # new winner, keep going
                    min_isl = out.nisl
                    opt_box = self.rms_box
                elif out.nisl == min_isl: # tie, keep going
                    pass
                else: # number of islands increased, stop
                    break
            else:
                # usually happens when box size is too small and "an
                # unphysical rms value was encountered", so stop here
                sf_logger.info('PyBDSF has failed with box {}.'.format(
                    self.rms_box))
                sf_logger.info('Stopping here.')
                break
            i += 1

        sf_logger.info('Found minimum of {} islands using box {}.'.format(
            min_isl, opt_box))
        # Final complete run with best parameters
        sf_logger.info('Final run...')
        self.rms_box = opt_box
        self.stop_at = None
        opt_out = self.find_sources()
        return opt_out
