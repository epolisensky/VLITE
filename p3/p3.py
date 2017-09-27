import os
import sys
from datetime import datetime
from errors import ConfigError
from database import createdb, dbclasses, filldb
from sourcefinding import runPyBDSF
from matching import radioXmatch
from skycatalog import catalogio
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
        s1 = (data['stages'])['save_to_database']
        s1qa = (data['stages'])['quality_checks']
        s2 = (data['stages'])['source_finding']
        s3 = (data['stages'])['source_association']
        s4 = (data['stages'])['catalog_matching']
        rootdir = (data['setup'])['rootdir']
        yr = (data['setup'])['year']
        mo = (data['setup'])['month']
        days = (data['setup'])['day']
        dbpath = (data['setup'])['dbname']
        owrite = (data['setup'])['overwrite']
        reproc = (data['setup'])['reprocess']
        catdb = (data['setup'])['catdb']
        catalogs = (data['setup'])['catalogs']
        params = data['pybdsf_params']
        res_tol = (data['matching'])['res_tol']
        catmatch = (data['matching'])['catalog_match']

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
    dirs = []
    for day in days:
        if type(day) != 'string':
            day = str(day)
        procdir = os.path.join(monthdir, day, 'Images/')
        if not os.path.isdir(procdir): # check full image path
            raise ConfigError('Directory does not exist: {}'.format(procdir))
        dirs.append(procdir)

    # Perform checks based on stages to process
    if s1 or s3:
        if dbpath is None:
            raise ConfigError('No database specified.')
        # Ensure reprocess & overwrite are boolean True/False or yes/no
        if isinstance(owrite, bool):
            pass
        else:
            raise ConfigError('Setup: overwrite must be True/False or yes/no.')
        if s3:
            if isinstance(reproc, bool):
                pass
            else:
                raise ConfigError('Setup: reprocess must be True/False '
                                  'or yes/no.')

    try:
        # Create & check path to sky catalogs database
        skycat = os.path.join(catalogio.catalogdir, catdb)
    except AttributeError:
        skycat = None
    if s4:
        if not os.path.isfile(skycat):
            raise ConfigError('Sky catalogs database does not exist: {}'.format(
            skycat))

        # Make sure requested catalogs exist
        catalog_opts = ['FIRST', 'GLEAM', 'NVSS', 'SUMSS', 'TGSS', 'WENSS']
        for cat in catalogs:
            if type(cat) != 'string':
                cat = str(cat)
            if cat not in catalog_opts:
                print('\nCurrently available catalogs: {}\n'.format(
                    catalog_opts))
                raise ConfigError('Catalog {} is not a valid option'.format(
                    cat))

        # Cross-match all or only new sources with sky catalogs
        catalog_match_opts = ['new', 'all']
        if catmatch not in catalog_match_opts:
            raise ConfigError('Matching: catalog_match must be all or new.')

    if s3 or s4:
        # Check resolution tolerance input
        try:
            res_tol = float(res_tol)
        except ValueError:
            raise ConfigError('Matching: res_tol must be a number.')
        except TypeError: # None
            pass

    stages = (s1, s1qa, s2, s3, s4)
      
    return stages, dirs, dbpath, owrite, skycat, catalogs, reproc, params, \
        res_tol, catmatch


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


def stage1(add, qa, dbname, impath):
    # STAGE 1 -- Initialize image object, add to table, 1st quality check
    # Add image to database?
    if add:
        # See if image is already in database
        exists = filldb.existCheck(dbname, impath)
        if exists is not None:
            if not reproc:
                imobj = None
            else:
                imobj = filldb.addImage(dbname, impath, exists)
        else:
            imobj = filldb.addImage(dbname, impath, exists)
    else:
        imobj = filldb.initImage(impath)

    # Run quality checks?
    if qa:
        # quality check
        pass
    else:
        pass
    
    return imobj


def stage2(add, qa, dbname, impath, params, imobj):
    # STAGE 2 -- Source finding + 2nd quality check
    print('\nExtracting sources from {}'.format(impath))
    # Initialize source finding image object
    bdsfim = runPyBDSF.BDSFImage(impath)
    for key, value in params.items():
        setattr(bdsfim, key, value)
    if params['mode'] == 'default':
        out = bdsfim.find_sources() # Run PyBDSF source finding
    else:
        out = bdsfim.minimize_islands()

    if out is not None:
        # Write PyBDSF(M) files to daily directory
        runPyBDSF.write_sources(out)
        imobj, sources = filldb.pipe_translate(imobj, out)
        if add:
            filldb.addSources(dbname, imobj, sources)
        if qa:
            # run more quality checks
            pass
        else:
            pass
    else:
        sources = None

    return imobj, sources


def stage3(dbname, imobj, sources, res_tol):
    # STAGE 3 - Source association
    radius = imobj.search_radius()
    center = (imobj.obs_ra, imobj.obs_dec)
    dassoc = radioXmatch.cone_search(dbname, ['AssocSource'], center, radius,
                                     same_res=True, beam=imobj.bmaj,
                                     res_tol=res_tol)
    #src = dassoc[dassoc.keys()[0]][0]
    #print src['ra']
    if not dassoc[dassoc.keys()[0]]: # table is empty, so add all sources
        assocsrcs = sources
        for asrc in assocsrcs:
            asrc.beam = imobj.bmaj
            asrc.num_detect = 1
            asrc.num_null = 0
        filldb.addAssoc(dbname, assocsrcs)
        for src in sources:
            src.assoc_id = src.src_id + 1
        filldb.updateAssocid(dbname, sources)
    else:
        print('\nMatching new sources to existing VLITE sources')
        


def stage4(imobj, sources, catalogs):
    # STAGE 4 - Catalog cross-matching
    matched_srcs, non_matched_srcs, matched_catsrcs = radioXmatch.main(
        catalogs=catalogs, database=False, objects=(imobj, sources))

    return matched_srcs, non_matched_srcs, matched_catsrcs

            
def process(stages, dirs, dbname, skycat, catalogs, reproc, params,
            res_tol, catmatch):   
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

            # STAGE 1 -- Initialize image object, add to table, QA
            imobj = stage1(stages[0], stages[1], dbname, impath)
            if imobj is None: # skip re-processing
                continue

            # STAGE 2 -- Source finding
            if stages[2]:
                imobj, sources = stage2(stages[0], stages[1], dbname, impath,
                                        params, imobj)
                if sources is None:
                    # Image failed to process
                    with open(pybdsfdir+'failed.txt', 'a') as f:
                        f.write(img+'\n')
                    if stages[0]:
                        filldb.pybdsf_fail(dbname, imobj)
                    continue

            # STAGE 3 -- Source association
            if stages[3]:
                stage3(dbname, imobj, sources, res_tol)


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

        #os.system('mv '+imgdir+'*pybdsf* '+pybdsfdir+'.')
        #os.system('mv '+imgdir+'*pybdsm* '+pybdsfdir+'.')




# Start the timer
start_time = datetime.now()

try:
    cf = sys.argv[1]
except IndexError:
    raise ConfigError('Please provide a configuration file.')
# Parse config file
stages, dirs, dbname, owrite, skycat, catalogs, reproc, params, rtol, \
    catmatch = cfgparse(cf)

if stages[0] or stages[3]:
    # Create/overwrite database if requested
    dbinit(dbname, owrite)

# Process images
process(stages, dirs, dbname, skycat, catalogs, reproc, params, rtol, catmatch)

print('\nTotal runtime: {}\n'.format(datetime.now() - start_time))
