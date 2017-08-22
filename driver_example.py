import os
import sys
import datetime
sys.path.insert(0, '/data3/erichards/codes/SourceFinding')
import runPyBDSF
import database
sys.path.insert(0, '/data3/erichards/codes/CatalogMatching')
import crossmatch


def print_info(start_time, outdir):
    """Print general information about the just
    completed run of the Post-Processing Pipeline."""
    # Count the number of images processed
    runPyBDSF.BDSFImage.image_count()
    # Count the number of catalogs written
    catnum = len([f for f in os.listdir(outdir) if f[-4:] == '.srl'])
    print ("\nWrote {} catalogs.\n".format(catnum))
    # Print the runtime
    print("\nTotal runtime: {}\n".format(datetime.datetime.now() - start_time)) 


# Start the timer
start_time = datetime.datetime.now()


# ========================================================
# User-defined parameters
# ========================================================
# Define path to images & database
rootdir = '/nfsshare/vpipe/processed/2017-08/01/'
database = '/data3/erichards/vlite/allvlite.sqlite'
# Define which and what order to check catalogs
catalogs = ['TGSS', 'NVSS', 'FIRST', 'SUMSS', 'WENSS']
# ========================================================

imgdir = os.path.join(rootdir, 'Images/')
pybdsfdir = os.path.join(rootdir, 'PyBDSF/')
if not os.path.isdir(pybdsfdir):
    os.system('mkdir '+pybdsfdir)
else:
    pass

# Create database if it doesn't already exist
#  - User will be prompted to verify creation
if not os.path.isfile(database):
    sourcedb.create_db(database)
else:
    pass

# Select only the images that end with 'IPln1.fits'
imglist = [f for f in os.listdir(imgdir) if f[-10:] == 'IPln1.fits']

# Begin loop through images
for img in imglist:
    impath = os.path.join(imgdir, img)
    # Define any desired PyBDSF parameters here
    bdsfim = runPyBDSF.BDSFImage(impath, trim_box=(400, 1300, 400, 1300),
                                 thresh='hard')
    out = bdsfim.find_sources() # Run PyBDSF source finding

    if out is not None:

        # Write PyBDSF(M) files to daily directory
        runPyBDSF.write_sources(out)
        os.system('mv '+imgdir+' *pybdsf* '+pybdsfdir+'.')
        os.system('mv '+imgdir+' *pybdsm* '+pybdsfdir+'.')

        '''Write sources to database; sources will be overwritten
        if this image has already been processed. The ImageTable
        and DetectedSources objects are returned to be used in
        cross-matching.'''
        imgtbl, sources = database.pybdsf_to_db(database, out)
        # Run catalog cross-matching
        ? = radioXmatch.main(catalogs=catalogs, database=False,
                             objects=(imgtbl, sources))
    else:
        # Image failed to process
        with open(pybdsfdir+'failed.txt', 'a') as f:
            f.write(img+'\n')


print_info(start_time, pybdsfdir)
