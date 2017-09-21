import os
import sys
from datetime import datetime
from errors import ConfigError
from database import createdb, dbclasses, filldb
from sourcefinding import runPyBDSF
from catalogmatching import radioXmatch
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


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
    print('\nTotal runtime: {}\n'.format(datetime.now() - start_time))


def cfgparse(cfgfile): 
    with open(cfgfile, 'r') as stream:
        data = load(stream, Loader=Loader)
        rootdir = (data['setup'])['rootdir']
        yr = (data['setup'])['year']
        mo = (data['setup'])['month']
        days = (data['setup'])['day']
        dbpath = os.path.join((data['setup'])['dbdir'],
                              (data['setup'])['dbname'])
        overwrite = (data['setup'])['overwrite']
        catalogs = (data['setup'])['catalogs']
        reprocess = (data['setup'])['reprocess']
        params = data['pybdsf_params']
        matchfor = (data['catalog_matching'])['do_for']

    # Perform checks on path to images
    if not os.path.isdir(rootdir): # check root directory
        raise ConfigError('Directory does not exist: {}'.format(rootdir))
    
    try:
        mo = format(int(mo), '02') # force 2-digits
    except ValueError: # will be caught in other exception
        pass

    yrmo = '{}-{}'.format(yr, mo)
    monthdir = os.path.join(rootdir, yrmo)
    if not os.path.isdir(monthdir): # check path with year-month
        raise ConfigError('Directory does not exist: {}'.format(monthdir))

    # If days = [], process every day in month directory
    if len(days) < 1:
        days = next(os.walk(monthdir))[1]

    # Define path to processing directories
    procdirs = []
    for day in days:
        if type(day) != 'string':
            day = str(day)
        procdir = os.path.join(monthdir, day, 'Images/')
        if not os.path.isdir(procdir): # check full image path
            raise ConfigError('Directory does not exist: {}'.format(procdir))
        procdirs.append(procdir)

    # Make sure requested catalogs exist
    catalog_opts = ['FIRST', 'GLEAM', 'NVSS', 'SUMSS', 'TGSS', 'WENSS']
    for cat in catalogs:
        if type(cat) != 'string':
            cat = str(cat)
        if cat not in catalog_opts:
            print('\nCurrently available catalogs: {}\n'.format(catalog_opts))
            raise ConfigError('Catalog {} is not a valid option'.format(cat))

    # Ensure reprocess & overwrite are boolean True/False or yes/no
    if isinstance(overwrite, bool):
        pass
    else:
        raise ConfigError('Setup: overwrite must be True/False or yes/no.')
    if isinstance(reprocess, bool):
        pass
    else:
        raise ConfigError('Setup: reprocess must be True/False or yes/no.')
    
    # Cross-match all or only new sources with sky catalogs
    catalog_match_opts = ['new', 'all']
    if matchfor not in catalog_match_opts:
        raise ConfigError('Catalog_matching: do_for must be all or new.')
    
    return procdirs, dbpath, overwrite, catalogs, reprocess, params, matchfor


def dbinit(dbname, overwrite):
    # Create database if it doesn't already exist
    #  - User will be prompted to verify creation
    if not os.path.isfile(dbname):
        makenew = raw_input('\nCreate new database {}? '.format(dbname))
        if makenew == 'y' or makenew == 'yes':
            safe = True
            createdb.create(dbname, safe)
        else:
            print('\nNo new database created.\n')
            raise ConfigError('Cannot access database {}'.format(dbname))
    else:
        if not overwrite:
            print('\nUsing existing database {}'.format(dbname))
        else:
            print('\nOverwriting existing database {}'.format(dbname))
            safe = False # this will prompt warning in create func
            createdb.create(dbname, safe)


def stage1(dbname, impath, exists):
    # STAGE 1 -- Initialize Image object & table + 1st quality check
    img = filldb.addImage(dbname, impath, exists)
    return img


def stage2(impath, params):
    # STAGE 2 - Source finding + 2nd quality check
    print('\nExtracting sources from {}'.format(impath))
    # Initialize source finding image object
    bdsfim = runPyBDSF.BDSFImage(impath)
    for key, value in params.items():
        setattr(bdsfim, key, value)
    if params['mode'] == 'default':
        out = bdsfim.find_sources() # Run PyBDSF source finding
    else:
        out = bdsfim.minimize_islands()

    return out


#def stage3():
    # STAGE 3 - Source association


def stage4(imobj, sources, catalogs):
    # STAGE 4 - Catalog cross-matching
    matched_srcs, non_matched_srcs, matched_catsrcs = radioXmatch.main(
        catalogs=catalogs, database=False, objects=(imobj, sources))

    return matched_srcs, non_matched_srcs, matched_catsrcs

            
def process(dirs, dbname, catalogs, reproc, params, matchfor):   
    # Begin loop through daily directories
    for imgdir in dirs:
        # Define/make directory for PyBDSF output
        daydir = os.path.abspath(os.path.join(imgdir, '..'))
        pybdsfdir = os.path.join(daydir, 'PyBDSF/')
        if not os.path.isdir(pybdsfdir):
            os.system('mkdir '+pybdsfdir)

        # Select only the images that end with 'IPln1.fits'
        imglist = [f for f in os.listdir(imgdir) if f.endswith('IPln1.fits')][:1]
        imglist.sort()

        # Begin loop through images
        for img in imglist:
            impath = os.path.join(imgdir, img)

            # See if image is already in database
            exists = filldb.existCheck(dbname, impath)
            if exists is not None:
                if not reproc:
                    break

            # 1.) 1st quality check + insertion into database
            imobj = stage1(dbname, impath, exists)

            # 2.) source finding
            out = stage2(impath, params)
            if out is not None:
                # Write PyBDSF(M) files to daily directory
                runPyBDSF.write_sources(out)
                os.system('mv '+imgdir+'*pybdsf* '+pybdsfdir+'.')
                os.system('mv '+imgdir+'*pybdsm* '+pybdsfdir+'.')

                '''Write sources to database; sources will be overwritten
                if this image has already been processed. The Image
                and DetectedSource objects are returned to be used in
                cross-matching. 2nd quality check is also performed here.'''
                imobj, sources = filldb.pipe_to_table(dbname, imobj, out)

                # 3.) Source association
                #stage3()
                
                # 4.) Catalog cross-matching
                '''
                if matchfor == 'new':
                    # extract only new sources from stage 3
                    pass
                else:
                    pass
                matched, not_matched, cat_matched = stage4(imobj, sources,
                                                           catalogs)
                # Write ds9 regions file of the matched catalog sources
                radioXmatch.write_region(cat_matched, impath,
                                         ext='_matches.reg')
                os.system('mv '+imgdir+'*.reg '+pybdsfdir+'.')
                # Update the database Source table with match info
                radioXmatch.update_database(dbname, imobj.filename, matched,
                                            not_matched)
                '''
            else:
                # Image failed to process
                filldb.pybdsf_fail(dbname, imobj)
                with open(pybdsfdir+'failed.txt', 'a') as f:
                    f.write(imobj+'\n')



# Start the timer
start_time = datetime.now()

cf = sys.argv[1]
# Parse config file
dirs, dbname, owrite, catalogs, reproc, params, matchfor = cfgparse(cf)

# Create/overwrite database if requested
dbinit(dbname, owrite)

# Process images
process(dirs, dbname, catalogs, reproc, params, matchfor)

print('\nTotal runtime: {}\n'.format(datetime.now() - start_time))
