"""runPyBDSF.py contains all the functionality to run a fits image
through the LOFAR source finding software PyBDSF 
(https://github.com/lofar-astron/PyBDSF). The class
BDSFImage creates an object whose attributes can be any number
of PyBDSF parameters and are specified through keyword arguments 
at initialization. Source finding is performed by calling the
object method find_sources(), which calls the PyBDSF function
process_image() and returns the output object. See sourcedb.py 
to see how to access results of the source finding as attributes 
of the bdsf output object.

Post-Processing Pipeline (P3) Stage 1"""


from astropy.io import fits
import numpy as np
import warnings
import bdsf
from database.dbclasses import Image


def write_sources(out):
    """Calls PyBDSF functions write_catalog() & export_image()
    to record results from source finding on a single image."""
    # Write the source list catalog, ascii and ds9 regions file
    out.write_catalog(format='ds9', catalog_type='srl', clobber=True)
    out.write_catalog(format='ascii', catalog_type='srl', clobber=True)

    # Write the residual image
    out.export_image(img_type='gaus_resid', clobber=True)
    # Write the model image
    # out.export_image(img_type='gaus_model', clobber=True)


class BDSFImage(Image):
    """Object to be manipulated and read into PyBDSF.
    Inherits all methods defined for Image class, but
    overrides initialization."""
 
    def __init__(self, image, box_incr=10, max_iter=10, **kwargs):
        """Initializes the Image subclass object. PyBDSF will 
        use any arguments it recognizes and ignore the rest."""
        self.filename = image
        self.box_incr = box_incr # used in minimize_islands
        self.max_iter = max_iter # used in minimize_islands
        self.quiet = True
        for key, value in kwargs.items():
            setattr(self, key, value)
        

    def fix_ctype3(self, data, header):
        """Change Obit-generated header keyword CTYPE3 from
        'SPECLNMF' to 'FREQ'."""
        try:
            if header['CTYPE3'] == 'SPECLNMF':
                header['CTYPE3'] = 'FREQ'
                self.write(data, header)
        except:
            try:
                data, header = self.read()
                if header['CTYPE3'] == 'SPECLNMF':
                    header['CTYPE3'] = 'FREQ'
                    self.write(data, header)
            except:
                print('\nERROR: Unable to update image header.\n')


    def get_attr(self):
        """Return all object attributes as a dictionary.
        This is fed into bdsf.process_image() as parameter
        arguments."""
        return self.__dict__


    def find_sources(self):
        """Run PyBDSF process_image() task using object attributes
        as parameter inputs. Returns None if PyBDSF fails."""
        opts = self.get_attr()
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', r'invalid value')
            try:
                out = bdsf.process_image(opts)
                return out
            except:
                return None  


    def minimize_islands(self):
        """Incrementally increase and decrease the bdsf.process_image()
        argument rms_box size until a minimum number of islands are 
        found. This technique works best on images with significant 
        artifacts where false detections are the biggest concern.
        This also takes a really long time since it is running
        bdsf.process_image numerous times, so use sparingly."""
        # Initial run
        self.stop_at = 'isl' # stop fitting at islands
        out = self.find_sources()
        if out is not None:
            box0 = out.rms_box # record initial box
            min_isl = out.nisl # initialize island number minimum
            print ("\nInitial box {} found {} islands.\n".format(box0,
                                                                  min_isl))
        else:
            # Only fails when user's box is too small
            box0 = self.rms_box            
            min_isl = 99999
            print ("\nPyBDSF has failed with box {}.".format(box0))
            print("Increasing rms_box...\n")

        opt_box = box0 # initialize optimal box
        box_size = box0[0]

        # Begin increasing box size loop
        i = 0
        while i < self.max_iter: # stop at max number of iterations
            box_size += self.box_incr
            box_step = int(box_size / 3.)
            self.rms_box = (box_size, box_step)
            print ("\nTrying box {}\n".format(self.rms_box))
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
                print ("\nPyBDSF has failed with box {}.".format(self.rms_box))
                print("Increasing rms_box...\n")
            i += 1

        # Begin decreasing box size loop
        i = 0
        box_size = box0[0] # reset the box size back to original
        while i < self.max_iter:
            box_size -= self.box_incr
            box_step = int(box_size / 3.)
            self.rms_box = (box_size, box_step)
            if box_size > 0: # stop if box size goes to 0 or below
                print ("\nTrying box {}\n".format(self.rms_box))
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
                print ("\nPyBDSF has failed with box {}.".format(self.rms_box))
                print ("Stopping here.\n")
                break
            i += 1

        print ("\nFound minimum of {} islands using box {}\n".format(min_isl,
                                                                    opt_box))
        # Final complete run with best parameters
        self.rms_box = opt_box
        self.stop_at = None
        opt_out = self.find_sources()
        return opt_out
