"""This is the main script for the VLITE Post-Processing 
Pipeline (P3). It is responsible for reading in the
configuration file, connecting to the the `PostgreSQL` database,
and calling the appropriate processing stages in the correct order.

"""
import os
import sys
import glob
import argparse
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from datetime import datetime
from errors import ConfigError
from database import createdb, dbclasses, dbio
from sourcefinding import runbdsf
from matching import radioxmatch
from skycatalog import catalogio, skycatdb
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


def print_run_stats(start_time):
    """Prints general information about the just completed run 
    of the Post-Processing Pipeline, including the number of
    images initialized and the total runtime.
    
    Parameters
    ----------
    start_time : datetime.datetime instance
        Time when the run was started.
    """
    print('--------------------------------------\n')
    print('Run statistics:')
    # Count the number of images processed
    nimgs = dbclasses.Image.image_count()
    print('\nProcessed {} images.'.format(nimgs))
    # Print the runtime
    runtime = datetime.now() - start_time
    print('\nTotal runtime: {}\n'.format(runtime))
    return nimgs, runtime


def cfgparse(cfgfile):
    """This function reads the `yaml` configuration file provided
    as a command line argument when `p3.py` is called.

    Parameters
    ----------
    cfgfile : str
        Name of the configuration file.

    Returns
    -------
    stages : tuple of 3 booleans
        Indicates which stages to run listed in the following order: 
        *source finding, source association, catalog matching*.
        Items can be ``True`` or ``False``.
    opts : tuple of 6 booleans
        Turns on and off options listed in the following order: 
        *save to database, quality checks, overwrite, reprocess, 
        redo match, update match*. Items can be ``True`` or ``False``.
    dirs : list of str
        List of strings specifying paths to daily image directories
        to be processed during the run.
    dbname : str
        Name of the `PostgreSQL` database.
    dbusr : str
        Username for the `PostgreSQL` database connection.
    catalogs : list of str
        List of sky survey tables to use when running catalog matching.
        They are cross-matched in the order they are listed.
    params : dict
        Specifies any non-default `PyBDSF` parameters to be used in source
        finding. 
    """    
    with open(cfgfile, 'r') as stream:
        data = load(stream, Loader=Loader)
        stages = data['stages']
        opts = data['options']
        setup = data['setup']
        sfparams = data['pybdsf_params']
        qaparams = data['image_qa_params']

    rootdir = setup['root directory']
    yr = setup['year']
    mo = setup['month']
    days = setup['day']
      
    if not any(stages.values()) and not opts['save to database']:
        raise ConfigError('Nothing to do -- change options: save to database: '
                          'to yes in the config file to write images to '
                          'database, or add a processing stage.')

    # Perform checks on path to images 
    try:
        mo = format(int(mo), '02') # force 2-digits
    except: # Value/TypeError will be caught in isdir exception below
        pass

    yrmo = '{}-{}'.format(yr, mo)
    try:
        monthdir = os.path.join(rootdir, yrmo)
    except AttributeError:
        raise ConfigError('Please provide a valid data root directory.')
    if not os.path.isdir(monthdir): # check path with year-month
        raise ConfigError('Directory does not exist: {}'.format(monthdir))

    # If days = [], process every day in month directory
    try:
        if len(days) < 1:
            days = next(os.walk(monthdir))[1]
    except TypeError:
        raise ConfigError('setup: day: must be a list (i.e. [{}])'.format(
            days))

    # Define path to processing directories
    dirs = []
    for day in days:
        if type(day) != 'string':
            day = str(day)
        procdir = os.path.join(monthdir, day, 'Images/')
        if not os.path.isdir(procdir): # check full image path
            raise ConfigError('Directory does not exist: {}'.format(procdir))
        dirs.append(procdir)

    # Make sure stage & option inputs are boolean
    for stage in stages.values():
        if isinstance(stage, bool):
            pass
        else:
            raise ConfigError('stage inputs must be True/False or yes/no.')
    for opt in opts.values():
        if isinstance(opt, bool):
            pass
        else:
            raise ConfigError('option inputs must be True/False or yes/no.')

    # Catch case when no DB is given
    if setup['database name'] is None or setup['database user'] is None:
        raise ConfigError('Please provide a database name/user.')

    if stages['catalog matching']:
        # Make sure requested catalogs exist
        catalog_opts = catalogio.catalog_list
        try:
            if len(setup['catalogs']) < 1:
                setup['catalogs'] = None
            for cat in setup['catalogs']:
                if type(cat) != 'string':
                    cat = str(cat)
                cat = cat.lower()
                if cat not in catalog_opts:
                    print('\nCurrently available catalogs: {}\n'.format(
                        catalog_opts))
                    raise ConfigError('Catalog {} is not a valid option'.format(
                        cat))
        except TypeError:
            raise ConfigError('Please provide a list of valid sky catalogs.')

    # Set default QA requirements if not specified
    if opts['quality checks']:
        if qaparams['min time on source (s)'] is None:
            qaparams['min time on source (s)'] = 60.
        else:
            try:
                qaparams['min time on source (s)'] = float(
                    qaparams['min time on source (s)'])
            except ValueError:
                raise ConfigError('min time on source must be a number.')
        if qaparams['max noise (mJy/beam)'] is None:
            qaparams['max noise (mJy/beam)'] = 1000.
        else:
            try:
                qaparams['max noise (mJy/beam)'] = float(
                    qaparams['max noise (mJy/beam)'])
            except ValueError:
                raise ConfigError('max noise must be a number.')
        if qaparams['max beam axis ratio'] is None:
            qaparams['max beam axis ratio'] = 4.
        else:
            try:
                qaparams['max beam axis ratio'] = float(
                    qaparams['max beam axis ratio'])
            except ValueError:
                raise ConfigError('max beam axis ratio must be a number.')
        if qaparams['min problem source separation (deg)'] is None:
            qaparams['min problem source separation (deg)'] = 20.
        else:
            try:
                qaparams['min problem source separation (deg)'] = float(
                    qaparams['min problem source separation (deg)'])
            except ValueError:
                raise ConfigError('min problem source separation must be a '
                                  'number.')
        if qaparams['max source metric'] is None:
            qaparams['max source metric'] = 10.
        else:
            try:
                qaparams['max source metric'] = float(
                    qaparams['max source metric'])
            except ValueError:
                raise ConfigError('max source metric must be a number.')
            
    return stages, opts, setup, sfparams, qaparams, dirs


def dbinit(dbname, user, overwrite, qaparams, safe_override=False):
    """Creates a `psycopg2` connection object to communicate
    with the `PostgreSQL` database. If no database with the
    provided name exists, the user is prompted to create a
    new one. This is done by first connecting to the *postgres* 
    database and then calling the SQL CREATE DATABASE command. 
    If the database is new or if overwriting the old one, all 
    necessary schemas, tables, and functions are additionally
    created through a call to `database.createdb.create()`.
    
    Parameters
    ----------
    dbname : str
        Name of the `PostgreSQL` database.
    user : str
        Username for the `PostgreSQL` database connection.
    overwrite : bool
        If ``True``, tables and data will be deleted and re-created,
        assuming there is a pre-existing database of name `dbname`.
        If ``False`, the existing database with `dbname` is used.
    safe_override : bool, optional
        If ``True``, this overrides the safe boolean. Implemented
        for testing purposes. Default value is ``False``.

    Returns
    -------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    """
    try:
        # DB exists
        conn = psycopg2.connect(host='localhost', database=dbname, user=user)
        print('\nConnected to database {}.'.format(dbname))
        if not overwrite:
            print('\nUsing existing database {}.'.format(dbname))
        else:
            print('\nOverwriting existing database tables.')
            if safe_override:
                createdb.create(conn, qaparams, safe=True)
            else:
                # This will prompt warning in create func
                createdb.create(conn, qaparams, safe=False)
        # Check for sky catalogs by verifying schema exists
        cur = conn.cursor()
        cur.execute('''SELECT EXISTS(SELECT 1 FROM pg_namespace 
            WHERE nspname = 'skycat');''')
        if not cur.fetchone()[0]:
            cur.close()
            print('\nNo sky catalogs found in database. '
                  'Creating them now...')
            skycatdb.create(conn)
        else:
            cur.close()
    except psycopg2.OperationalError:
        # DB does not yet exist
        if safe_override:
            makenew = 'yes'
        else:
            makenew = raw_input('\nCreate new database {}? '.format(dbname))
        if makenew == 'y' or makenew == 'yes':
            conn = psycopg2.connect(host='localhost', database='postgres',
                                    user=user)
            conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
            cur = conn.cursor()
            cur.execute('CREATE DATABASE ' + dbname)
            cur.close()
            conn.close()
            conn = psycopg2.connect(host='localhost', database=dbname,
                                    user=user)
            print('\nConnected to new database {}.'.format(dbname))
            createdb.create(conn, qaparams, safe=True)
            print('\nCreating new sky catalog tables...')
            skycatdb.create(conn)
        else:
            print('\nNo new database created.\n')
            raise ConfigError('Cannot access database {}'.format(dbname))

    return conn


def iminit(conn, impath, save, qa, qaparams, reproc, stages):
    """Initializes `database.dbclasses.Image()` object
    and adds row to database image table if the image is not
    already in the table or if the image is being reprocessed.
    This function represents the first stage in the 
    Post-Processing Pipeline.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    impath : str
        Directory path to the FITS image file.
    save : bool
        If ``True``, the image info is written and saved
        to the database image table. If ``False``, an Image
        class object is initialized without writing to the database.
    qa : bool
        Turns on and off quality assurance. If ``True``, quality
        checks are run on the image data.
    reproc : bool
        Turns on and off reprocessing. If ``True``, the image is
        re-initialized as an Image class object and updated in the
        database image table (if save = ``True``) even if the image
        already has an entry in the database. If ``False`` and the
        image is in the database, no further processing is done and
        the code moves on to the next image in the list.

    Returns
    -------
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info, or ``None`` if the image
        has already been processed and will not be reprocessed.
    """
    # STAGE 1 -- Initialize image object, add to table, 1st quality check

    # Initialize Image object & set attributes from header
    imobj = dbio.init_image(impath)

    # Is the image in the database?
    status = dbio.status_check(conn, impath)

    # Run image quality checks only if it is a new image
    if qa:
        imobj.image_qa(qaparams)
    else:
        pass

    # Only adding image to DB - add or update w/o deleting sources
    if not any(stages.values()):
        if status is None:
            # image not in DB
            imobj = dbio.add_image(conn, imobj, status) # branch 2
        else:
            if reproc:
                # already processed, but re-doing -- sources NOT deleted
                imobj = dbio.add_image(conn, imobj, status) # branch 4
            else:
                # already processed & not re-doing
                print('\nImage {} already in database. Moving on...'.format(
                    impath))
                imobj = None # branch 3
    else: # Running at least one stage
        if status is None:
            if stages['source finding']:
                if save:
                    # image not in DB; planning to source find & write to DB
                    imobj = dbio.add_image(conn, imobj, status)
                else:
                    # just initialize if not writing to DB
                    print('\nInitializing {}'.format(impath))
            else:
                # image not in DB, but not running source finding --> quit
                print('\nERROR: Image {} not yet processed. Source finding '
                      'must be run before other stages.'.format(impath))
                imobj = None # branch 5
        else:
            if stages['source finding']:
                if reproc:
                    if save:
                        # image in DB; re-doing SF, so old sources are removed
                        imobj = dbio.add_image(conn, imobj, status,
                                               delete=True)
                    else:
                        # just initialize if not writing to DB
                        print('\nInitializing {}'.format(impath))
                else:
                    # image already in DB & not re-processing
                    print('\nImage {} already processed. Moving on...'.format(
                        impath))
                    imobj = None # branch 9
            else:
                # not running SF -- stage must be > 1
                if status[1] > 1:
                    print('\nInitializing {}'.format(impath))
                    imobj.id = status[0]
                    imobj.stage = status[1]
                    imobj.radius = status[2]
                else:
                    print('\nERROR: Image {} does not have sources extracted '
                          'yet. Source finding must be run before other '
                          'stages.'.format(impath))
                    imobj = None # branch 7
    
    return imobj


def vlite_unique(conn, src, image_id, radius):
    existing = dbio.check_vlite_unique(conn, src.id)
    if not existing: # returned empty list
        src.detected = True
        # Add source to vlite_unique table
        dbio.add_vlite_unique(conn, src, image_id)
        # Find previous images where VU source is in FOV
        src.detected = False
        prev_images = radioxmatch.check_previous(conn, src, radius)
        for previd in prev_images:
            # Ignore current image
            if previd[0] != image_id:
                dbio.add_vlite_unique(conn, src, previd[0])
    else: # VU source already in table
        existing_imgids = [row[1] for row in existing]
        if image_id not in existing_imgids:
            # Add if new image
            src.detected = True
            dbio.add_vlite_unique(conn, src, image_id)
        else: # entry exists
            entry = [row for row in existing if row[1] == image_id]
            if not entry[0][2]:
                # Update entry to detected = True
                src.detected = True
                dbio.add_vlite_unique(conn, src, image_id, update=True)

    return src


def srcfind(conn, imobj, sfparams, save, qa, qaparams):
    """Runs `PyBDSF` source finding, writes out results,
    and inserts sources into database detected_source and
    detected_island tables, if the save to database option is on.
    This function represents the second stage of the
    Post-Processing Pipeline.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info.
    params : dict
        Specifies any non-default `PyBDSF` parameters to be 
        used in source finding.
    save : bool
        If ``True``, the extracted sources are written and saved
        to the database detected_island and detected_source tables. The
        `PyBDSF` files are always written out.
    qa : bool
        Turns on and off quality assurance. If ``True``, quality
        checks are run on the source finding results.

    Returns
    -------
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    sources : list
        List of database.dbclasses.DetectedSource objects.
        Attributes of each object are set from the `PyBDSF`
        output object.
    """
    # STAGE 2 -- Source finding + 2nd quality check
    print('\nExtracting sources from {}'.format(imobj.filename))
    # Initialize source finding image object
    bdsfim = runbdsf.BDSFImage(imobj.filename, **sfparams)
    # Run PyBDSF source finding
    if sfparams['mode'] == 'minimize_islands':
        out = bdsfim.minimize_islands()
    else:
        out = bdsfim.find_sources()

    # Set new radius attribute for Image object & update stage
    imobj.radius = bdsfim.radius
    imobj.stage = 2
   
    if out is not None:
        # Write PyBDSF files to daily directory
        runbdsf.write_sources(out)
        imobj, sources = dbclasses.translate(imobj, out)
        if qa:
            imobj.source_qa(sources, qaparams)
    else:
        # PyBDSF failed to process
        sources = None
        imobj.error_id = 5

    if save:
        dbio.add_sources(conn, imobj, sources)
        # Compute beam corrected fluxes & write to corrected_flux table
        print('\nCorrecting all flux measurements for primary beam response')
        for src in sources:
            src.correct_flux(imobj)
            dbio.add_corrected(conn, src)

    return imobj, sources


def srcassoc(conn, imobj, sources, save):
    """Associates sources extracted from current image with
    previously detected VLITE sources stored in the
    assoc_source database table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    sources : list
        List of database.dbclasses.DetectedSource objects.
        Attributes of each object are set from the `PyBDSF`
        output object.
    save : bool
        If ``True``, the assoc_source table is updated with the
        association results and the assoc_id is updated in the
        detected_source table. An entry is also added to the
        vlite_unique table if a source with 0 catalog matches
        is pulled from the asssoc_source table.
    
    Returns
    -------
    detected_unmatched : list of DetectedSource objects
        Sources extracted from the image which could not be successfully
        associated with previously detected sources. These are added
        to the assoc_table as new detections.
    imobj : database.dbclasses.Image instance
        Initialized Image object with updated stage attribute.
    """
    # STAGE 3 -- Source association
    detected_matched, detected_unmatched, assoc_matched, assoc_unmatched \
        = radioxmatch.associate(conn, sources, imobj, imobj.radius)
    if save:
        # Update assoc_id col for matched detected source
        if detected_matched:
            dbio.update_detected_associd(conn, detected_matched)
        # Add new (unmatched) detected sources to assoc_source table
        if detected_unmatched:
            # Updates assoc_id attribute
            detected_unmatched = dbio.add_assoc(conn, detected_unmatched)
        # Update matched assoc_source entries
        if assoc_matched:
            dbio.update_matched_assoc(conn, assoc_matched)
        # Check for VLITE unique (VU) sources that weren't detected in image
        for asrc in assoc_unmatched:
            if asrc.nmatches == 0:
                asrc.detected = False
                # Add source to vlite_unique table
                dbio.add_vlite_unique(conn, asrc, imobj.id)
        # Check for VU sources that were detected in image
        for asrc in assoc_matched:
            if asrc.nmatches == 0:
                asrc.detected = True
                dbio.add_vlite_unique(conn, asrc, imobj.id)

        # Update stage in image table
        imobj.stage = 3
        dbio.update_stage(conn, imobj)

    return assoc_matched, detected_unmatched, imobj


def catmatch(conn, imobj, sources, catalogs, save):
    """Cross-matches a list of DetectedSource objects
    to all specified sky catalogs.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    sources : list
        List of database.dbclasses.DetectedSource objects.
    catalogs : list of str
        Names of the sky survey catalogs to use.
    save : bool
        If ``True``, match results are recorded in the
        catalog_match table and the assoc_source table is updated.
        VLITE unique sources with no sky catalog match are
        also inserted into the vlite_unique table.
    """
    # STAGE 4 -- Sky catalog cross-matching
    # Filter out catalogs which have already been checked for this image
    new_catalogs = dbio.update_checked_catalogs(conn, imobj.id, catalogs)
    for catalog in new_catalogs:
        # Filter out sources which already have results for this catalog
        todo_sources = []
        for src in sources:
            already_matched = dbio.check_catalog_match(conn, src.id, catalog)
            if not already_matched:
                todo_sources.append(src)
        try:
            todo_sources, catalog_matched = radioxmatch.catalogmatch(
                conn, todo_sources, catalog, imobj, imobj.radius)
        except TypeError: # no sky catalog sources extracted
            continue
        if save:
            # Add results to catalog_match table
            dbio.add_catalog_match(conn, catalog_matched)
    if save:
        # Update assoc_source.nmatches
        dbio.update_assoc_nmatches(conn, todo_sources)
        # Check for new VLITE unique (VU) sources from this image
        for src in todo_sources:
            if src.nmatches == 0:
                src  = vlite_unique(conn, src, imobj.id, imobj.radius)
                            
        # Update stage
        imobj.stage = 4
        dbio.update_stage(conn, imobj)

    return

            
def process(conn, stages, opts, dirs, catalogs, sfparams, qaparams):
    """
    This function handles the logic and transitions between
    processing stages. All individual processing functions
    are called from here. Input parameters come from output
    of `cfgparse`.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    stages : tuple of 3 booleans
        Indicates which stages to run listed in the following order: 
        *source finding, source association, catalog matching*.
        Items can be ``True`` or ``False``.
    opts : tuple of 6 booleans
        Turns on and off options listed in the following order: 
        *save to database, quality checks, overwrite, reprocess,
        redo match, update match*. Items can be ``True`` or ``False``.
    dirs : list of str
        List of strings specifying paths to daily image directories
        to be processed during the run.
    catalogs : list of str
        List of sky survey tables to use when running catalog matching.
    sfparams : dict
        Specifies any non-default `PyBDSF` parameters to be used in source
        finding. 
    """
    sf = stages['source finding']
    sa = stages['source association']
    cm = stages['catalog matching']
    save = opts['save to database']
    qa = opts['quality checks']
    reproc = opts['reprocess']
    rematch = opts['redo match']
    updatematch = opts['update match']
    
    # Begin loop through daily directories
    for imgdir in dirs:
        # Define/make directory for PyBDSF output
        daydir = os.path.abspath(os.path.join(imgdir, '..'))
        pybdsfdir = os.path.join(daydir, 'PyBDSF/')
        if not os.path.isdir(pybdsfdir):
            os.system('mkdir '+pybdsfdir)

        # Select only the images that end with 'IPln1.fits'
        imglist = [f for f in os.listdir(imgdir) if f.endswith('IPln1.fits')]
        imglist.sort()
        #imglist = imglist[:1]

        # Begin loop through images
        for img in imglist:
            impath = os.path.join(imgdir, img)

            # STAGE 1 -- Initialize image object
            imobj = iminit(conn, impath, save, qa, qaparams, reproc, stages)
            # Move on to next image if imobj is None
            if imobj is None:
                continue
            
            # STAGE 2 -- Source finding
            if sf:
                imobj, sources = srcfind(
                    conn, imobj, sfparams, save, qa, qaparams)
                os.system('rm '+imgdir+'*.crop.fits')
                if sources is None:
                    # Image failed to process
                    with open(pybdsfdir+'failed.txt', 'a') as f:
                        f.write(img+'\n')
                    continue
                else:
                    os.system('mv '+imgdir+'*pybdsf.log '+pybdsfdir+'.')
                    if glob.glob(imgdir+'*pybdsm.srl'):
                        os.system('mv '+imgdir+'*pybdsm* '+pybdsfdir+'.')
                # STAGE 3 -- Source association
                if sa:
                    matched_assoc, new_sources, imobj = srcassoc(
                        conn, imobj, sources, save)
                    # STAGE 4 -- Sky survey catalog cross-matching
                    if cm: # sf + sa + cm - branch 12, 15
                        # Cross-match new sources only
                        catmatch(conn, imobj, new_sources, catalogs, save)
                        if glob.glob(imgdir+'*matches.reg'):
                            os.system(
                                'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                        print('\n======================================\n')
                        print('Completed source finding, association, '
                              'and sky catalog cross-matching.')
                        print('\n======================================\n')
                    else: # sf + sa - branch 11, 14
                        print('\n======================================\n')
                        print('Completed source finding and association.')
                        print('\n======================================\n')
                        continue
                else:
                    if cm: # sf + cm - branch 10, 13
                        catmatch(conn, imobj, sources, catalogs, False)
                        if glob.glob(imgdir+'*matches.reg'):
                            os.system(
                                'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                        print('\n======================================\n')
                        print('Completed source finding and sky catalog '
                              'cross-matching to the extracted sources.')
                        print('\n======================================\n')
                    else: # sf only - branch 6, 8
                        print('\n======================================\n')
                        print('Completed source finding.')
                        print('\n======================================\n')
                        continue
            else: # no sf
                if sa:
                    # Get sources from detected_source table
                    sources = dbio.get_image_sources(conn, imobj.id)
                    # Already caught case of no sf but stage < 2 in iminit
                    if imobj.stage == 2: # no sa has been run yet
                        matched_assoc, new_sources, imobj = srcassoc(
                            conn, imobj, sources, save)
                        if cm: # sa + cm - branch 20
                            # Cross-match new sources only
                            catmatch(conn, imobj, new_sources, catalogs, save)
                            if glob.glob(imgdir+'*matches.reg'):
                                os.system(
                                    'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                            print('\n======================================\n')
                            print('Completed source association and sky '
                                  'catalog cross-matching to the newly '
                                  'detected sources.')
                            print('\n======================================\n')
                        else: # sa only - branch 19
                            print('\n======================================\n')
                            print('Completed associating detected sources with '
                                  'the existing catalog.')
                            print('\n======================================\n')
                    else: # stage > 2
                        print("\n{}'s sources have already been associated "
                              "with the existing catalog.".format(impath))
                        if cm: # cm only - branch 21
                            assoc_sources = dbio.get_associated(conn, sources)
                            if rematch:
                                assoc_sources = dbio.delete_matches(
                                    conn, assoc_sources, imobj.id)
                            else:
                                if not updatematch:
                                    assoc_sources = [src for src in \
                                                     assoc_sources if \
                                                     src.nmatches is None \
                                                     or src.nmatches == 0]
                                else: pass
                            catmatch(conn, imobj, assoc_sources, catalogs, save)
                            if glob.glob(imgdir+'*matches.reg'):
                                os.system(
                                    'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                            print('\n======================================\n')
                            print('Completed sky catalog cross-matching.')
                            print('\n======================================\n')
                        else: continue # branch 18
                else:
                    if cm: # cm only - branch 17
                        pass
                    else: # branches 2, 4 continued
                        continue
                    if imobj.stage > 2:
                        # Get sources from detected_source table
                        sources = dbio.get_image_sources(conn, imobj.id)
                        assoc_sources = dbio.get_associated(conn, sources)
                        if rematch:
                            assoc_sources = dbio.delete_matches(
                                conn, assoc_sources, imobj.id)
                        else:
                            if not updatematch:
                                assoc_sources = [src for src in assoc_sources \
                                                 if src.nmatches is None or \
                                                 src.nmatches == 0]
                            else: pass
                        catmatch(conn, imobj, assoc_sources, catalogs, save)
                        if glob.glob(imgdir+'*matches.reg'):
                            os.system(
                                'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                        print('\n======================================\n')
                        print('Completed sky catalog cross-matching.')
                        print('\n======================================\n')
                    else: # branch 16
                        print('\nERROR: Source association for image {} must '
                              'be run before catalog cross-matching. '
                              'Alternatively, run source finding first and '
                              'cross-match extracted sources with sky survey '
                              'catalogs.'.format(impath))
                        continue

    return

            
def main():
    """One function to rule them all."""
    parser = argparse.ArgumentParser(
        description='Run the VLITE database Post-Processing Pipline (P3)')
    parser.add_argument('config_file', help='the YAML configuration file')
    parser.add_argument('--ignore_prompt', action='store_true',
                        help='ignore prompt to verify database '
                        'removal/creation')
    parser.add_argument('--remove_catalog_matches', action='store_true',
                        help='remove matching results for the specified '
                        'sky survey catalog(s)')
    parser.add_argument('--remove_source', action='store_true',
                        help='removes the specified source(s) from the '
                        'database assoc_source table')
    parser.add_argument('--add_catalog', action='store_true',
                        help='adds any new sky survey catalogs to a table in '
                        'the database "skycat" schema')
    args = parser.parse_args()
    
    # Start the timer
    start_time = datetime.now()

    stages, opts, setup, sfparams, qaparams, dirs = cfgparse(args.config_file)

    # Find existing/create/overwrite database
    conn = dbinit(setup['database name'], setup['database user'],
                  opts['overwrite'], qaparams, safe_override=args.ignore_prompt)

    if args.remove_catalog_matches:
        catalogs = raw_input('\nFor which catalogs would you like to remove '
                             'matching results? (List catalogs separated by '
                             'a comma.)\n')
        cat_list = [cat.lower() for cat in catalogs.split(', ')]
        print('\nRemoving matching results for {}...'.format(catalogs))
        dbio.remove_catalog(conn, cat_list)
        # Find all assoc_sources whose nmatches dropped to 0
        vu_assoc_sources = dbio.get_new_vu(conn)
        if vu_assoc_sources is not None:
            for vu_asrc in vu_assoc_sources:
                # Get image_ids, radii for the new VU source
                vu_image_list = dbio.get_vu_image(conn, vu_asrc.id)
                for vu_image in vu_image_list:
                    src = vlite_unique(conn, vu_asrc, vu_image[0], vu_image[1])
        conn.close()
        sys.exit(0)

    if args.remove_source:
        inp = raw_input('\nPlease enter the id number(s) (separated by '
                        'commas) of the source(s) you wish to remove '
                        'from the database assoc_source table or provide '
                        'a path to a text file with the id numbers:\n')
        try:
            asid_list = [int(asid) for asid in inp.strip('[]').split(', ')]
        except ValueError:
            with open(inp, 'r') as f:
                text = f.read()
            asid_list = [int(asid) for asid in text.strip().split('\n')]
        print('\nRemoving row(s) {} from the assoc_source table...'.format(
            asid_list))
        dbio.remove_sources(conn, tuple(asid_list))
        conn.close()
        sys.exit(0)

    if args.add_catalog:
        skycatdb.create(conn)
        conn.close()
        sys.exit(0)

    # Process images
    process(conn, stages, opts, dirs, setup['catalogs'], sfparams, qaparams)

    nimages, exec_time = print_run_stats(start_time)

    # Record run configuration parameters
    dbio.record_config(conn, args.config_file, start_time, exec_time,
                       nimages, stages, opts, setup, sfparams, qaparams)

    conn.close()


if __name__ == '__main__':
    main()
