"""This is the main script for the VLITE Database
Pipeline (vdp). It is responsible for reading in the
configuration file, connecting to the the `PostgreSQL` database,
and calling the processing stages.

"""
import os
import sys
import glob
import re
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


__version__ = 1.8


def print_run_stats(start_time):
    """Prints general information about the just completed run 
    of the Post-Processing Pipeline, including the number of
    images initialized and the total execution time.
    
    Parameters
    ----------
    start_time : datetime.datetime instance
        Time when the run was started.

    Returns
    -------
    nimgs : integer
        Number of image files initialized as Image objects.
    runtime : datetime.timedelta instance
        Total execution time of the just completed pipeline run.
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
    """This function reads the `YAML` configuration file provided
    as a command line argument when `vdp.py` is called and parses
    each section of the file into dictionaries.

    Parameters
    ----------
    cfgfile : str
        Name of the configuration file.

    Returns
    -------
    stages : dict
        Keys are the processing stages (source finding, source assocation,
        and catalog matrching) and values are boolean ``True`` or ``False``.
    opts : dict
        Keys are the processing options (save to database, quality checks,
        overwrite, reprocess, redo match, and update match) and values are
        boolean ``True`` or ``False``.
    setup : dict
        Keys are the setup parameters (root directory, year, month, day,
        database name, database user, and catalogs) and values are the
        user-supplied inputs.
    sfparams : dict
        Keys are the source finding (mode and scale) and `PyBDSF` parameters
        and values are the user inputs.
    qaparams : dict
        Keys are the image quality parameters (min time on source (s),
        max noise (mJy/beam), max beam axis ratio, min problem source
        separation (deg), and max source metric) and values are the user
        inputs. Defaults are defined if none are specified.        
    dirs : list
        List of strings specifying paths to daily image directories
        to be processed during the run.
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
    days = sorted(setup['day'])

    # Raise error if the configuration says to do nothing
    if not any(stages.values()) and not opts['save to database']:
        raise ConfigError('Nothing to do -- change options: save to database: '
                          'to "yes" in the configuration file to write images '
                          'to the database or add a processing stage.')

    # Perform checks on path to images
    if yr is None or mo is None:
        procdir = os.path.join(rootdir, 'Images/')
        if not os.path.isdir(procdir):
            raise ConfigError('Directory does not exist: {}'.format(procdir))
        dirs = [procdir]
    else:
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
                days = sorted(next(os.walk(monthdir))[1])
        except TypeError:
            raise ConfigError('setup: day: must be a list (i.e. [{}])'.format(
                days))

        # Define path to processing directories
        dirs = []
        for day in days:
            try:
                day = format(int(day), '02')
            except ValueError:
                continue
            procdir = os.path.join(monthdir, day, 'Images/')
            # Check full image path
            if not os.path.isdir(procdir):
                print('\nSkipping non-existent directory {}'.format(procdir))
                continue
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

    # Catch case when no database is given
    if setup['database name'] is None or setup['database user'] is None:
        raise ConfigError('Please provide a database name/user.')
    # Force database name to all lowercase
    setup['database name'] = setup['database name'].lower()

    # Check list of sky catalogs
    if stages['catalog matching']:
        catalog_opts = catalogio.catalog_dict.keys() # all available catalogs
        # If catalogs = [], use all of them
        if len(setup['catalogs']) < 1:
            setup['catalogs'] = catalog_opts
        else:
            # Make sure requested catalogs exist
            try:
                for cat in setup['catalogs']:
                    if type(cat) != str:
                        cat = str(cat)
                    cat = cat.lower()
                    if cat not in catalog_opts:
                        print('\nCurrently available catalogs: {}\n'.format(
                            catalog_opts))
                        raise ConfigError('Catalog {} is not a valid option'.
                                          format(cat))
            except TypeError:
                raise ConfigError('Please provide a list of valid sky '
                                  'catalogs.')

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
    new one. The "skycat" schema which holds all the sky survey
    catalogs in tables is created at this stage if it does not
    already exist. The user will be prompted to verify deletion
    of all current tables if the database exists and the *overwrite*
    option in the configuration file is ``True``. All necessary
    tables, functions, and triggers are created through a call to 
    `database.createdb.create()`.
    
    Parameters
    ----------
    dbname : str
        Name of the `PostgreSQL` database.
    user : str
        Username for the `PostgreSQL` database connection.
    overwrite : bool
        If ``True``, tables and data will be deleted and re-created,
        assuming there is a pre-existing database of name "dbname".
        If ``False``, the existing database with "dbname" is used.
    qaparams : dict
        The dictionary of image quality cuts specified by the user
        in the configuration file. The values are written to the
        database **error** table.
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
                # This will prompt warning in create function
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


def vlite_unique(conn, src, image_id, radius):
    """This function adds a VLITE unique (VU) source, or a source
    with no sky catalog matches, to the **vlite_unique** table. 
    After adding a new VU source, the **image** table is queried
    to find all previously processed images in which the VU source
    could have been detected based on field-of-view (FOV) and spatial
    resolution. The **vlite_unique** table, therefore, keeps a record
    of every image in which a VU source was in the FOV and whether it
    was detected.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    src : database.dbclasses.DetectedSource instance
        The VU source from the **assoc_source** table.
    image_id : int
        Id number of the image in which the VU source
        was detected.
    radius : float
        Radius (in degrees) of the image's FOV. Used in
        querying the **image** table.

    Returns
    -------
    src : database.dbclasses.DetectedSource instance
        The same VU src with updated `detected` attribute.
    """
    # Start by checking if the VU source is already in the VU table
    existing = dbio.check_vlite_unique(conn, src.id)
    if not existing: # returned empty list
        # Newly detected VU source
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
    else:
        # VU source already in table, check if new image
        existing_imgids = [row[1] for row in existing]
        if image_id not in existing_imgids:
            # Add if new image
            src.detected = True
            dbio.add_vlite_unique(conn, src, image_id)
        else:
            # previously undetected VU source, now detected
            # (i.e. after re-doing source finding)
            entry = [row for row in existing if row[1] == image_id]
            if not entry[0][2]: # if not detected
                # Update entry to detected = True
                src.detected = True
                dbio.add_vlite_unique(conn, src, image_id, update=True)

    return src


def iminit(conn, imobj, save, qa, qaparams, reproc, stages):
    """This function handles ingestion of the image metadata
    into the database **image** table and represents the first
    stage in the Post-Processing Pipeline.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from the header info.
    save : bool
        If ``True``, the image info is written and saved
        to the database **image** table.
    qa : bool
        If ``True``, quality checks are run on the image data.
    qaparams : dict
        Dictionary of image quality requirements and their
        user-specified values from the configuration file.
    reproc : bool
        If ``True``, any existing entry for this image in the
        database **image** table is updated. If ``False`` and
        there is an existing entry, no further processing is
        done and the code moves on to the next image in the list.
    stages : dict
        Dictionary specifying which processing stages are to be run.

    Returns
    -------
    imobj : database.dbclasses.Image instance
        Image object with `id` attribute updated after insertion into
        the database **image** table, or ``None`` if the image
        has already been processed and will not be reprocessed.
    """
    # STAGE 1 -- Add image to table, 1st quality check
    print
    print('**********************')
    print('STAGE 1: READING IMAGE')
    print('**********************')

    # Run image quality checks
    if qa:
        imobj.image_qa(qaparams)
    else:
        pass

    # Is the image in the database?
    status = dbio.status_check(conn, imobj.filename)

    # Only adding image to DB - add or update w/o deleting sources
    if not any(stages.values()):
        if status is None:
            # image not in DB, so add it
            imobj = dbio.add_image(conn, imobj, status) # branch 2
        else:
            if reproc:
                # already processed, but re-doing -- sources NOT deleted
                imobj = dbio.add_image(conn, imobj, status) # branch 4
            else:
                # already processed & not re-doing
                print('\nImage already in database. Moving on...')
                imobj = None # branch 3
    else: # Running at least one stage
        if status is None:
            if stages['source finding']:
                if save:
                    # image not in DB; planning to source find & write to DB
                    imobj = dbio.add_image(conn, imobj, status)
                else:
                    # just initialize if not writing to DB
                    print('\nInitializing image.')
            else:
                # image not in DB, but not running source finding --> quit
                print('\nERROR: Image {} not yet processed. Source finding '
                      'must be run before other stages.'.format(imobj.filename))
                imobj = None # branch 5
        else:
            if stages['source finding']:
                if status[1] == 1:
                    # image in DB, but no SF yet
                    if save:
                        imobj = dbio.add_image(conn, imobj, status)
                    else:
                        print('\nInitializing image.')
                else:
                    if reproc:
                        if save:
                            # image in DB, delete old SF results
                            imobj = dbio.add_image(conn, imobj, status,
                                                   delete=True)
                        else:
                            # just initialize if not writing to DB
                            print('\nInitializing image.')
                    else:
                        # image has SF results & not re-processing
                        print('\nImage already processed. Moving on...')
                        imobj = None # branch 9
            else:
                # not running SF -- stage must be > 1
                if status[1] > 1:
                    print('\nInitializing image.')
                    imobj.id = status[0]
                    imobj.stage = status[1]
                    imobj.radius = status[2]
                else:
                    print('\nERROR: Image {} does not have sources extracted '
                          'yet. Source finding must be run before other '
                          'stages.'.format(imobj.filename))
                    imobj = None # branch 7

    # Uncomment below when ready to stop flagged images (after testing)
    #if imobj.error_id is not None:
    # For now, just stop if header keywords are missing
    if imobj is not None and imobj.error_id == 1:
        imobj = None
    
    return imobj


def srcfind(conn, imobj, sfparams, save, qa, qaparams):
    """Runs `PyBDSF` source finding, writes out results,
    and inserts source fit parameters into database 
    **detected_source**, **detected_island**, and 
    **corrected_flux** tables, if the *save to database*
    option is set to ``True``. This function represents
    the second stage of the Post-Processing Pipeline.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info.
    sfparams : dict
        Specifies any non-default `PyBDSF` parameters to be 
        used in source finding.
    save : bool
        If ``True``, the source fit parameters are written and saved
        to the database **detected_island*, **detected_source**, and
        **corrected_flux** tables. The `PyBDSF` files are always 
        written out.
    qa : bool
        If ``True``, quality checks are run on the source finding
        results.
    qaparams : dict
        User-specified requirements from the configuration file for
        the source finding quality checks.

    Returns
    -------
    imobj : database.dbclasses.Image instance
        Initialized Image object with updated attributes
        from the source finding results.
    sources : list
        List of database.dbclasses.DetectedSource objects.
        Attributes of each object are set from the `PyBDSF`
        output object.
    """
    # STAGE 2 -- Source finding + 2nd quality check
    print
    print('**********************')
    print('STAGE 2: SOURCE FINDNG')
    print('**********************')
    
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
        # Translate `PyBDSF` output to DetectedSource objects
        imobj, sources = dbclasses.translate(imobj, out)
        # Run quality checks, part 2
        if qa:
            imobj.source_qa(sources, qaparams)
    else:
        # PyBDSF failed to process
        sources = None
        imobj.error_id = 6

    if save:
        # Add source fit parameters to database tables
        dbio.add_sources(conn, imobj, sources)
        if sources is not None:
            # Compute beam corrected fluxes & write to corrected_flux table
            print('Correcting all flux measurements for primary beam '
                  'response.')
            for src in sources:
                src.correct_flux(imobj)
                dbio.add_corrected(conn, src)

    return imobj, sources


def srcassoc(conn, imobj, sources, save):
    """Associates through positional cross-matching sources
    extracted from the current image with previously detected
    VLITE sources stored in the **assoc_source** database table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    sources : list
        List of database.dbclasses.DetectedSource objects.
        Attributes of each object are from the `PyBDSF` fit results.
    save : bool
        If ``True``, the **assoc_source** table is updated with the
        association results and the 'assoc_id' is updated in the
        **detected_source** table. If ``False``, no results are
        saved to the database.
    
    Returns
    -------
    detected_unmatched : list
        List of new VLITE detected sources. 
    imobj : database.dbclasses.Image instance
        Initialized Image object with updated `stage` attribute.
    """
    # STAGE 3 -- Source association
    print
    print('***************************')
    print('STAGE 3: SOURCE ASSOCIATION')
    print('***************************')

    if save:
        match_in_db = True
    else:
        match_in_db = False

    # Associate current sources with existing VLITE catalog 
    detected_matched, detected_unmatched, assoc_matched, assoc_unmatched \
        = radioxmatch.associate(conn, sources, imobj, imobj.radius, match_in_db)
    if save:
        # Update assoc_id col for matched detected sources
        if detected_matched:
            dbio.update_detected_associd(conn, detected_matched)
        # Add new (unmatched) detected sources to assoc_source table
        if detected_unmatched:
            # Updates assoc_id attribute
            detected_unmatched = dbio.add_assoc(conn, detected_unmatched)
        # Update matched assoc_source positions
        if assoc_matched:
            dbio.update_matched_assoc(conn, assoc_matched)
        # Check for VLITE unique (VU) sources that weren't detected in image
        for asrc in assoc_unmatched:
            if asrc.nmatches == 0:
                asrc.detected = False
                # Update vlite_unique table with image/source non-detection
                dbio.add_vlite_unique(conn, asrc, imobj.id)
        # Check for VU sources that were detected in image
        for asrc in assoc_matched:
            if asrc.nmatches == 0:
                asrc.detected = True
                dbio.add_vlite_unique(conn, asrc, imobj.id)

        # Update stage in image table
        imobj.stage = 3
        dbio.update_stage(conn, imobj)

    return detected_unmatched, imobj


def catmatch(conn, imobj, sources, catalogs, save):
    """Performs positional cross-matching of VLITE 
    detected sources to other radio sky survey catalogs.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info.
    sources : list
        VLITE detected sources to be matched to sky catalog sources.
    catalogs : list
        Names of the sky survey catalogs to use.
    save : bool
        If ``True``, match results are recorded in the
        **catalog_match** table and the **assoc_source** table is
        updated. VLITE unique sources with no sky catalog match are
        inserted into the **vlite_unique** table.
    """
    # STAGE 4 -- Sky catalog cross-matching
    print
    print('*********************************')
    print('STAGE 4: MATCHING TO SKY CATALOGS')
    print('*********************************')
    
    catalogs = [catalog.lower() for catalog in catalogs]
    # Filter catalogs by resolution
    filtered_catalogs = radioxmatch.filter_catalogs(conn, catalogs, imobj.bmin)
    # Remove catalogs that have already been checked for this image
    if save:
        new_catalogs = dbio.update_checked_catalogs(
            conn, imobj.id, filtered_catalogs)
    else:
        new_catalogs = catalogs
    if not new_catalogs:
        print('\nAll specified catalogs with appropriate resolution '
              'have already been checked for matches.')
        if save:
            imobj.stage = 4
            dbio.update_stage(conn, imobj)
        return

    print('\nUsing the following catalogs for cross-matching: {}'.format(
        new_catalogs))

    if not sources:
        print('\nNo new VLITE sources to match.')
        if save:
            imobj.stage = 4
            dbio.update_stage(conn, imobj)
        return

    # Cross-match VLITE sources with each catalog
    if not save:
        match_in_db = False
    else:
        match_in_db = True
    for catalog in new_catalogs:
        try:
            sources, catalog_matched = radioxmatch.catalogmatch(
                conn, sources, catalog, imobj, imobj.radius, match_in_db)
        except TypeError:
            # No sky catalog sources extracted, move on to next catalog
            continue
        if save:
            # Add results to catalog_match table
            dbio.add_catalog_match(conn, catalog_matched)
    if save:
        # Update assoc_source nmatches
        dbio.update_assoc_nmatches(conn, sources)
        # Check for new VLITE unique (VU) sources from this image
        for src in sources:
            if src.nmatches == 0:
                src = vlite_unique(conn, src, imobj.id, imobj.radius)
                            
        # Update stage
        imobj.stage = 4
        dbio.update_stage(conn, imobj)                

    return


def process(conn, stages, opts, dirs, files, catalogs, sfparams, qaparams):
    """This function handles the logic and transitions between
    processing stages. All individual processing functions
    are called from here. Input parameters come from output
    of `cfgparse`.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    stages : dict
        Keys are the processing stages (source finding, source assocation,
        and catalog matrching) and values are boolean ``True`` or ``False``.
    opts : dict
        Keys are the processing options (save to database, quality checks,
        overwrite, reprocess, redo match, and update match) and values are
        boolean ``True`` or ``False``.
    dirs : list
        List of strings specifying paths to daily image directories
        to be processed during the run.
    files : list
        List of files to process in each daily directory.
    catalogs : list
        Names of radio sky survey catalogs to use when running
        catalog matching.
    sfparams : dict
        Specifies any non-default `PyBDSF` parameters to be used in source
        finding.
    qaparams : dict
        User-specified quality requirements or default values defined
        and set in `cfgparse`.
    """
    # Define booleans from stages & opts dictionaries
    sf = stages['source finding']
    sa = stages['source association']
    cm = stages['catalog matching']
    save = opts['save to database']
    qa = opts['quality checks']
    reproc = opts['reprocess']
    rematch = opts['redo match']
    updatematch = opts['update match']
    
    # Begin loop through daily directories
    i = 0
    for imgdir in dirs:
        # Define/make directory for PyBDSF output
        daydir = os.path.abspath(os.path.join(imgdir, '..'))
        pybdsfdir = os.path.join(daydir, 'PyBDSF/')
        if not os.path.isdir(pybdsfdir):
            os.system('mkdir '+pybdsfdir)

        if not files[0]:
            # Select all images that end with 'IPln1.fits'
            imglist = [f for f in os.listdir(imgdir) if \
                       f.endswith('IPln1.fits')]
        else:
            imglist = [f for f in files[i]]
        i += 1

        # Loop through images to initialize
        imobjlist = []
        for img in imglist:
            impath = os.path.join(imgdir, img)
            # Initialize Image object & set attributes from header
            imobjlist.append(dbio.init_image(impath))

        # Sort imobjlist by mjdtime
        imobjlist.sort(key=lambda x: x.mjdtime)

        # Begin loop through time-sorted images
        for imobj in imobjlist:
            print('_' * (len(imobj.filename) + 10))
            print('\nStarting {}.'.format(imobj.filename))
            # STAGE 1 -- Add image to database
            imobj = iminit(conn, imobj, save, qa, qaparams, reproc, stages)
            # Move on to next image if imobj is None
            if imobj is None:
                continue
            
            # STAGE 2 -- Source finding
            if sf:
                imobj, sources = srcfind(
                    conn, imobj, sfparams, save, qa, qaparams)
                # Remove temp fits file
                os.system('rm '+imgdir+'*.crop.fits')
                # Do this so PyBDSF warning messages get printed to log
                # when re-directing out put to text file
                os.system('grep "WARNING" '+imgdir+'*.log')
                if sources is None:
                    # Image failed to process
                    with open(pybdsfdir+'failed.txt', 'a') as f:
                        f.write(img+'\n')
                    continue
                else:
                    # Move PyBDSF output files to PyBDSF directory
                    os.system('mv '+imgdir+'*pybdsf.log '+pybdsfdir+'.')
                    if glob.glob(imgdir+'*pybdsm*'):
                        os.system('mv '+imgdir+'*pybdsm* '+pybdsfdir+'.')
                # STAGE 3 -- Source association
                if sa:
                    new_sources, imobj = srcassoc(conn, imobj, sources, save)
                    # STAGE 4 -- Sky survey catalog cross-matching
                    if cm: # sf + sa + cm - branch 12, 15
                        # Cross-match new sources only
                        catmatch(conn, imobj, new_sources, catalogs, save)
                        if glob.glob(imgdir+'*matches.reg'):
                            os.system(
                                'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                        print('\n====================================='
                              '==========================================\n')
                        print('Completed source finding, association, '
                              'and sky catalog cross-matching on image')
                        print('{}.'.format(imobj.filename))
                        print('\n====================================='
                              '==========================================\n')
                    else: # sf + sa - branch 11, 14
                        print('\n====================================='
                              '==========================================\n')
                        print('Completed source finding and association on '
                              'image')
                        print('{}.'.format(imobj.filename))
                        print('\n====================================='
                              '==========================================\n')
                        continue
                else:
                    if cm: # sf + cm - branch 10, 13
                        catmatch(conn, imobj, sources, catalogs, False)
                        if glob.glob(imgdir+'*matches.reg'):
                            os.system(
                                'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                        print('\n====================================='
                              '==========================================\n')
                        print('Completed source finding and sky catalog '
                              'cross-matching to the extracted')
                        print('sources from image')
                        print('{}.'.format(imobj.filename))
                        print('\n====================================='
                              '==========================================\n')
                    else: # sf only - branch 6, 8
                        print('\n====================================='
                              '==========================================\n')
                        print('Completed source finding on image')
                        print('{}.'.format(imobj.filename))
                        print('\n====================================='
                              '==========================================\n')
                        continue
            else: # no sf
                if sa:
                    # Get sources from detected_source table
                    sources = dbio.get_image_sources(conn, imobj.id)
                    # Already caught case of no sf but stage < 2 in iminit
                    if imobj.stage == 2: # no sa has been run yet
                        new_sources, imobj = srcassoc(conn, imobj,
                                                      sources, save)
                        if cm: # sa + cm - branch 20
                            # Cross-match new sources only
                            catmatch(conn, imobj, new_sources, catalogs, save)
                            if glob.glob(imgdir+'*matches.reg'):
                                os.system(
                                    'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                            print('\n====================================='
                              '==========================================\n')
                            print('Completed source association and sky '
                                  'catalog cross-matching to the newly')
                            print('detected sources from image')
                            print('{}.'.format(imobj.filename))
                            print('\n====================================='
                              '==========================================\n')
                        else: # sa only - branch 19
                            print('\n====================================='
                              '==========================================\n')
                            print('Completed source association for image')
                            print('{}.'.format(imobj.filename))
                            print('\n====================================='
                              '==========================================\n')
                    else: # stage > 2
                        print("\nNOTE: {}'s".format(imobj.filename))
                        print('sources have already been associated with the '
                              'existing VLITE catalog.')
                        if cm: # cm only - branch 21
                            assoc_sources = dbio.get_associated(conn, sources)
                            if rematch:
                                # Delete & redo matching
                                assoc_sources = dbio.delete_matches(
                                    conn, assoc_sources, imobj.id)
                            else:
                                if not updatematch:
                                    # Cross-match new/un-matched sources only
                                    assoc_sources = [src for src in \
                                                     assoc_sources if \
                                                     src.nmatches is None \
                                                     or src.nmatches == 0]
                                else:
                                    # Use all sources if updating
                                    pass
                            catmatch(conn, imobj, assoc_sources, catalogs, save)
                            if glob.glob(imgdir+'*matches.reg'):
                                os.system(
                                    'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                            print('\n====================================='
                              '==========================================\n')
                            print('Completed sky catalog cross-matching for '
                                  'image')
                            print('{}.'.format(imobj.filename))
                            print('\n====================================='
                              '==========================================\n')
                        else: continue # branch 18
                else:
                    if cm: # cm only - branch 17
                        pass
                    else: # branches 2, 4 continued
                        continue
                    if imobj.stage > 2:
                        # Get detected, then assoc sources
                        sources = dbio.get_image_sources(conn, imobj.id)
                        assoc_sources = dbio.get_associated(conn, sources)
                        if rematch:
                            # Delete & redo matching
                            assoc_sources = dbio.delete_matches(
                                conn, assoc_sources, imobj.id)
                        else:
                            if not updatematch:
                                # Cross-match new/un-matched sources only
                                assoc_sources = [src for src in assoc_sources \
                                                 if src.nmatches is None or \
                                                 src.nmatches == 0]
                            else:
                                # Use all sources if updating
                                pass
                        catmatch(conn, imobj, assoc_sources, catalogs, save)
                        if glob.glob(imgdir+'*matches.reg'):
                            os.system(
                                'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                        print('\n====================================='
                              '==========================================\n')
                        print('Completed sky catalog cross-matching for image')
                        print('{}.'.format(imobj.filename))
                        print('\n====================================='
                              '==========================================\n')
                    else: # branch 16
                        print('\n====================================='
                              '==========================================\n')
                        print('\nERROR: Source association must be run '
                              'before catalog cross-matching for image')
                        print('{}.\n'.format(imobj.filename))
                        print('Alternatively, run source finding with catalog '
                              'matching to cross-match detected')
                        print('VLITE sources with other radio sky surveys.')
                        print('\n====================================='
                              '==========================================\n')
                        continue

    return


def main():
    """One function to rule them all."""
    # Set required & optional command line arguments
    parser = argparse.ArgumentParser(
        description='Run the VLITE Database Pipline (vdp)')
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
    parser.add_argument('--remove_image', action='store_true',
                        help='removes the specified image(s) and associated '
                        'results from the database entirely')
    parser.add_argument('--manually_add_match', action='store_true',
                        help='manually add catalog matching results for '
                        'VLITE source(s) after follow-up')
    parser.add_argument('--add_catalog', action='store_true',
                        help='adds any new sky survey catalogs to a table in '
                        'the database "skycat" schema')
    args = parser.parse_args()
    
    # Start the timer
    start_time = datetime.now()

    # Parse run configuration file
    stages, opts, setup, sfparams, qaparams, dirs = cfgparse(args.config_file)

    # Find existing/create/overwrite database
    if any([args.remove_catalog_matches, args.remove_source,
            args.remove_image, args.manually_add_match, args.add_catalog]):
        opts['overwrite'] = False
    conn = dbinit(setup['database name'], setup['database user'],
                  opts['overwrite'], qaparams, safe_override=args.ignore_prompt)

    # Option to remove matching results for sky catalogs
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

    # Option to remove sources from the **assoc_source** database table
    if args.remove_source:
        inp = raw_input('\nPlease enter the id number(s) (i.e. 1, 2, 3) '
                        'of the source(s) you wish to remove from the '
                        'database assoc_source table, or provide a text '
                        'file with one id number per line:\n')
        try:
            asid_list = [int(asid) for asid in inp.strip('[]').split(',')]
        except ValueError:
            with open(inp, 'r') as f:
                text = f.read()
            asid_list = [int(asid) for asid in text.strip().split('\n')]
        print('\nRemoving row(s) {} from the assoc_source table...'.format(
            asid_list))
        dbio.remove_sources(conn, tuple(asid_list))
        conn.close()
        sys.exit(0)

    # Option to remove images
    if args.remove_image:
        inp = raw_input('\nPlease enter the image(s) filename(s) starting '
                        'at least with the year-month directory (i.e. '
                        '2018-01/15/Images/10GHz.Mrk110.IPln1.fits), or '
                        'provide a text file with one filename per line:\n')
        try:
            images = [re.findall('([0-9]{4}-\S+)', img)[0] for img \
                      in inp.split(',')]
        except IndexError:
            with open(inp, 'r') as f:
                text = f.read()
            images = [re.findall('([0-9]{4}-\S+)', img)[0] for img \
                      in text.strip().split('\n')]
        print('\nPreparing to remove image(s) {} from the database.'.format(
            images))
        confirm = raw_input('\nAre you sure? ')
        if confirm == 'y' or confirm == 'yes':
            print('Deleting image(s) from the database...')
            dbio.remove_images(conn, images)
        else:
            print('Doing nothing...')
        conn.close()
        sys.exit(0)

    # Option to manually add catalog matching results
    if args.manually_add_match:
        inp = raw_input('\nPlease enter the source assoc_source id, the '
                        'name of the catalog, and, optionally, the id of the '
                        'matched catalog source and the angular separation in '
                        'arcseconds, in that order one per line. '
                        'Hit "q" when you are done. You may alternatively '
                        'provide a similarly formatted text file with one '
                        'catalog match per line:\n')
        cmatches = []
        while inp != 'q':
            try:
                int(inp[0])
                cmatches.append(inp)
                inp = raw_input()
            except IndexError:
                inp = raw_input()
            except ValueError:
                with open(inp, 'r') as f:
                    text = f.read()
                    cmatches = [line for line in text.strip().split('\n')]
                    break
            if inp == 'q':
                break
        assoc_ids = []
        catalogs = []
        catsrc_ids = []
        separations = []
        for cmatch in cmatches:
            cm = cmatch.split(', ')
            assoc_ids.append(int(cm[0]))
            catalog = cm[1].lower()
            if catalog not in catalogio.catalog_dict.keys():
                raise ConfigError('{} is not a valid catalog.'.format(cm[1]))
            else:
                catalogs.append(catalog)
            try:
                catsrc_ids.append(int(cm[2]))
                separations.append(float(cm[3]))
            except IndexError:
                catsrc_ids.append(-1)
                separations.append(-1)
        cmrows = zip(catalogs, catsrc_ids, assoc_ids, separations)
        print('\nAdding new catalog matching results for assoc_ids {}...'.
              format(assoc_ids))
        dbio.update_assoc_nmatches(conn, assoc_ids)
        dbio.add_catalog_match(conn, cmrows)
        conn.close()
        sys.exit(0)
    
    # Option to add a new sky survey catalog to the database "skycat" schema
    if args.add_catalog:
        skycatdb.create(conn)
        conn.close()
        sys.exit(0)
    
    # Process images
    process(conn, stages, opts, dirs, setup['files'],
            setup['catalogs'], sfparams, qaparams)

    nimages, exec_time = print_run_stats(start_time)

    # Record run configuration parameters
    dbio.record_config(conn, args.config_file, start_time, exec_time,
                       nimages, stages, opts, setup, sfparams, qaparams)

    conn.close()


if __name__ == '__main__':
    main()
