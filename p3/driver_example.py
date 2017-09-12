import os
import sys
import datetime
from database import createdb, dbclasses, filldb
from sourcefinding import runPyBDSF
from catalogmatching import radioXmatch


def print_run_stats(start_time, outdir):
    """Print general information about the just
    completed run of the Post-Processing Pipeline."""
    print('\n======================================\n')
    print('Run statistics:')
    # Count the number of images processed
    dbclasses.Image.image_count()
    # Count the number of catalogs written
    catnum = len([f for f in os.listdir(outdir) if f[-4:] == '.srl'])
    print ('\nWrote {} catalogs.'.format(catnum))
    # Print the runtime
    print('\nTotal runtime: {}\n'.format(datetime.datetime.now() - start_time)) 


# Start the timer
start_time = datetime.datetime.now()

rootpath = '/nfsshare/vpipe/processed/'
dbdir = '/data3/vpipe/'

# ========================================================
# User-defined parameters
# ========================================================
rootdir = os.path.join(rootpath, sys.argv[1])
dbname = os.path.join(dbdir, sys.argv[2])
catalogs = [c for c in sys.argv[3].split(', ')]
print('\nUsing {} for cross-matching.'.format(catalogs))
# ========================================================

# Catch potential user input errors
imgdir = os.path.join(rootdir, 'Images/')
if not os.path.isdir(imgdir):
    print('\nERROR: Image directory does not exist.')
    sys.exit(0)
else:
    pass

# Make sure given catalogs exist
catalog_opts = ['FIRST', 'GLEAM', 'NVSS', 'SUMSS', 'TGSS', 'WENSS']
for cat in catalogs:
    if cat not in catalog_opts:
        print('\nERROR: Catalog {} is not a valid option.'.format(cat))
        sys.exit(0)
    else:
        pass

# Create database if it doesn't already exist
#  - User will be prompted to verify creation
if not os.path.isfile(dbname):
    print('\nDatabase {} does not yet exist.'.format(dbname))
    print('It is safe to make a new one.')
    safe = True
    createdb.create(dbname, safe)
else:
    pass

# Define/make directory for PyBDSF output
pybdsfdir = os.path.join(rootdir, 'PyBDSF/')
if not os.path.isdir(pybdsfdir):
    os.system('mkdir '+pybdsfdir)
else:
    pass

# Select only the images that end with 'IPln1.fits'
imglist = [f for f in os.listdir(imgdir) if f[-10:] == 'IPln1.fits']

# Begin loop through images
for img in imglist:
    impath = os.path.join(imgdir, img)

    # STAGE 1 -- Initialize Image object & table
    # STAGE 2 -- Run preliminary quality checks
    # *********************************************
    img = filldb.addImage(dbname, impath)


    # STAGE 3 - Source finding
    # *******************************************
    # Define any desired PyBDSF parameters here
    bdsfim = runPyBDSF.BDSFImage(impath, thresh='hard')
    # trim_box=(400, 1300, 400, 1300)
    out = bdsfim.find_sources() # Run PyBDSF source finding

    if out is not None:
        # Write PyBDSF(M) files to daily directory
        runPyBDSF.write_sources(out)
        os.system('mv '+imgdir+'*pybdsf* '+pybdsfdir+'.')
        os.system('mv '+imgdir+'*pybdsm* '+pybdsfdir+'.')

        '''Write sources to database; sources will be overwritten
        if this image has already been processed. The Image
        and DetectedSources objects are returned to be used in
        cross-matching.'''
        img, sources = filldb.pipe_to_table(dbname, img, out)

        # STAGE 4 - Final image quality check
        # *******************************************

        # STAGE 5 - Source association
        # *******************************************
        

        # STAGE 6 - Catalog cross-matching
        # *******************************************
        # Run catalog cross-matching
        matched_srcs, non_matched_srcs, matched_catsrcs = radioXmatch.main(
            catalogs=catalogs, database=False, objects=(img, sources))
        # Write ds9 regions file of the matched catalog sources
        radioXmatch.write_region(matched_catsrcs, impath, ext='_matches.reg')
        os.system('mv '+imgdir+'*.reg '+pybdsfdir+'.')
        # Update the database Source table with match info
        radioXmatch.update_database(dbname, img.filename,
                                    matched_srcs, non_matched_srcs)

    else:
        # Image failed to process
        filldb.pybdsf_fail(dbname, img)
        with open(pybdsfdir+'failed.txt', 'a') as f:
            f.write(img+'\n')


print_run_stats(start_time, pybdsfdir)
