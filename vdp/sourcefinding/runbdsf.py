"""This module contains all the functionality to run a FITS image
through the LOFAR source finding software PyBDSF 
(https://github.com/lofar-astron/PyBDSF). The class
``BDSFImage()`` creates an object whose attributes can be any number
of PyBDSF parameters. Source finding is performed by calling the
object method ``find_sources()``, which calls the PyBDSF function
``process_image()`` and returns the output object.

"""
import warnings
from datetime import datetime
import numpy as np
import bdsf
from database.dbclasses import Image
from timeout import timeout


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
        self.scale = 1.0
        # Set our own default rms_box parameter
        self.set_rms_box()
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
        # Pre-process the image for PyBDSF
        self.crop()


    def crop(self):
        """Crops the FITS image into a circle so that sources
        are extracted from a region the same size and shape
        as when running a cone search query on the database.
        This ensures that the exact same areas on the sky are
        being considered when selecting sources for matching.

        """
        data, hdr = Image.read(self)
        Image.header_attrs(self, hdr)
        # Fix header keywords for PyBDSF
        if hdr['CTYPE3'] == 'SPECLNMF':
            hdr['CTYPE3'] = 'FREQ'
        try:
            hdr['BMAJ']
        except KeyError:
            hdr['BMAJ'] = self.bmaj / 3600. # deg
            hdr['BMIN'] = self.bmin / 3600.
            hdr['BPA'] = self.bpa
        n = len(data[0,0,:,:])
        a, b = hdr['CRPIX1'], hdr['CRPIX2']
        y, x = np.ogrid[-a:n-a, -b:n-b]
        r = (hdr['NAXIS2'] / 2.) * self.scale
        mask = x*x + y*y >= r*r
        self.radius = round(r * hdr['CDELT2'], 2) # in deg
        data[0,0,mask] = np.nan
        self.filename = self.filename[:-4]+'crop.fits'
        Image.write(self, data, hdr, owrite=True)


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
        print('\nExtracting sources...')
        start = datetime.now()
        opts = self.get_attr()
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', r'invalid value')
            try:
                out = bdsf.process_image(opts)
            except:
                out = None

        try:
            print(' -- found {} sources in {:.2f} seconds\n'.format(
                out.nsrc, (datetime.now() - start).total_seconds()))
        except AttributeError:
            try:
                print(' -- found {} islands in {:.2f} seconds\n'.format(
                    out.nisl, (datetime.now() - start).total_seconds()))
            except AttributeError:
                print(' -- PyBDSF failed to process image.\n')

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
        print('\nStarting minimize_islands with box {}...'.format(
            self.rms_box))
        self.stop_at = 'isl' # stop fitting at islands
        out = self.find_sources()
        if out is not None:
            box0 = out.rms_box # record initial box
            min_isl = out.nisl # initialize island number minimum
            #print ('Initial box {} found {} islands.\n'.format(box0, min_isl))
        else:
            # Only fails when user's box is too small
            box0 = self.rms_box            
            min_isl = 99999
            print ('\nPyBDSF has failed with box {}.'.format(box0))
            print('Increasing rms_box...\n')

        opt_box = box0 # initialize optimal box
        box_size = box0[0]

        # Begin increasing box size loop
        i = 0
        while i < self.max_iter: # stop at max number of iterations
            box_size += self.box_incr
            box_step = int(box_size / 3.)
            self.rms_box = (box_size, box_step)
            print ('Trying box {}...'.format(self.rms_box))
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
                print ('\nPyBDSF has failed with box {}.'.format(self.rms_box))
                print('Increasing rms_box...\n')
            i += 1

        # Begin decreasing box size loop
        i = 0
        box_size = box0[0] # reset the box size back to original
        while i < self.max_iter:
            box_size -= self.box_incr
            box_step = int(box_size / 3.)
            self.rms_box = (box_size, box_step)
            if box_size > 0: # stop if box size goes to 0 or below
                print ('Trying box {}...'.format(self.rms_box))
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
                print ('\nPyBDSF has failed with box {}.'.format(self.rms_box))
                print ('Stopping here.\n')
                break
            i += 1

        print('Found minimum of {} islands using box {}.'.format(min_isl,
                                                                    opt_box))
        # Final complete run with best parameters
        print('\nFinal run...')
        self.rms_box = opt_box
        self.stop_at = None
        opt_out = self.find_sources()
        return opt_out
