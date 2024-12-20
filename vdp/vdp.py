"""This is the main script for the VLITE Database
Pipeline (vdp). It is responsible for reading in the
configuration file, connecting to the the PostgreSQL database,
and calling the processing stages.

"""
import os
import sys
import glob
import re
import argparse
import logging
import psycopg2
import numpy as np
import healpy as hp
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from datetime import datetime
from break_handler import BreakHandler
from database import createdb, dbclasses, dbio
from errors import ConfigError
from sourcefinding import runbdsf,beam_tools
from matching import radioxmatch
from radiocatalogs import catalogio, radcatdb
from math import sqrt
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader
from astropy.io import fits

__version__ = '3.5'


# Create logger
logger = logging.getLogger('vdp')
logger.setLevel(logging.DEBUG)


def loggerinit(logfile=None, quiet=False):
    """Initializes handlers for logging to both the 
    console and a text file.

    Parameters
    ----------
    logfile : str, optional
        Name of the log file. If ``None``, messages will
        only be printed to the console. Default is ``None``.
    quiet : bool, optional
        If ``True``, no messages will be printed to the
        console. Default is ``False``.
    """
    # Create boolean lists to determine if any handlers of either type
    is_stream = []
    is_file = []
    if len(logger.handlers) > 0:
        for handler in logger.handlers:
            is_stream.append(type(handler) is logging.StreamHandler)
            is_file.append(type(handler) is logging.FileHandler)

    if quiet:
        pass
    elif any(is_stream):
        # a console handler already exists
        pass
    else:
        # create console handler with a slightly higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create formatter and add it to the handler
        ch_formatter = logging.Formatter('%(message)s')
        ch.setFormatter(ch_formatter)
        # add the handler to the logger
        logger.addHandler(ch)

    if logfile is None:
        pass
    elif any(is_file):
        # a file handler already exists
        pass
    else:
        # create file handler which logs even debug messages
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        # create formatter and add it to the handler
        fh_formatter = logging.Formatter(
            '%(asctime)s %(name)-28s: %(levelname)-6s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        fh.setFormatter(fh_formatter)
        # add the handler to the logger
        logger.addHandler(fh)


def print_run_stats(start_time):
    """Prints general information about the just completed run 
    of the pipeline, including the number of
    images processed and the total execution time.

    Parameters
    ----------
    start_time : ``datetime.datetime`` instance
        Time when the run was started.

    Returns
    -------
    nimgs : int
        Number of image files initialized as Image objects.
    runtime : ``datetime.timedelta`` instance
        Total execution time of the just completed pipeline run.
    """
    logger.info('--------------------------------------')
    logger.info('Run statistics:')
    # Count the number of images processed
    nimgs = dbclasses.Image.image_count()
    logger.info('Processed {} images.'.format(nimgs))
    # Print the runtime
    runtime = datetime.now() - start_time
    logger.info('Total runtime: {}'.format(runtime))
    return nimgs, runtime


def cfgparse(cfgfile):
    """This function reads the YAML configuration file provided
    as a command line argument when vdp.py is called and parses
    each section of the file into dictionaries.

    Parameters
    ----------
    cfgfile : str
        Name of the configuration file.

    Returns
    -------
    stages : dict
        Keys are the processing stages (source finding, source assocation,
        and catalog matching) and values are boolean ``True`` or ``False``.
    opts : dict
        Keys are the processing options (save to database, quality checks,
        overwrite, reprocess, redo match, and update match) and values are
        boolean ``True`` or ``False``.
    setup : dict
        Keys are the setup parameters (root directory, year, month, day,
        files, database name, database user, catalogs, & smear time) and values
        are the user-supplied inputs.
    sfparams : dict
        Keys are the required source finding parameters *mode*, *scale*,
        *borderpad* and other optional PyBDSF parameters and values are 
        the user inputs.
    qaparams : dict
        Keys are the image quality parameters (min time on source (s),
        max noise (mJy/beam), max beam axis ratio, min problem source
        separation (deg), max source metric, etc) and values are the user
        inputs. Defaults are defined if none are specified.        
    dirs : list
        List of strings specifying paths to daily image directories
        to be processed during the run.
    catflag : boolean
        Set True if user supplied catalogs for matching. Will match to all catalogs
        regardless of resolution.
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
    imgdir = setup['image directory']

    # Raise error if the configuration says to do nothing
    if not any(stages.values()) and not opts['save to database']:
        raise ConfigError('Nothing to do -- change options: save to database: '
                          'to "yes" in the configuration file to write images '
                          'to the database or add a processing stage.')

    # Perform checks on path to images
    if yr is None or mo is None:
        if not imgdir:
            imgdir='Images/'
        procdir=rootdir+imgdir
        print(procdir)
        if not os.path.isdir(procdir):
            raise ConfigError('Directory does not exist: {}'.format(procdir))
        if not days:
            dirs = [procdir]
        else:
            # Define path to processing directories
            dirs = []
            for day in days:
                procdir = rootdir+imgdir+day
                # Check full image path
                if not os.path.isdir(procdir):
                    print('\nSkipping non-existent directory {}'.format(procdir))
                    continue
                else:
                    dirs.append(procdir)
    else:
        try:
            mo = format(int(mo), '02')  # force 2-digits
        except:  # Value/TypeError will be caught in isdir exception below
            pass

        yrmo = '{}-{}'.format(yr, mo)
        try:
            monthdir = os.path.join(rootdir, yrmo)
        except AttributeError:
            raise ConfigError('Please provide a valid data root directory.')
        if not os.path.isdir(monthdir):  # check path with year-month
            raise ConfigError('Directory does not exist: {}'.format(monthdir))

        # If days = [], process every day in month directory
        try:
            if len(days) < 1:
                tmp = sorted(next(os.walk(monthdir))[1])
                for i in tmp:
                    if len(i) == 2:
                        days.append(i)  # assume only want the 2-digit day dirs
        except TypeError:
            raise ConfigError('setup: day: must be a list (i.e. [{}])'.format(
                days))

        # Set image directory to Images/ if left blank
        if not imgdir:
            imgdir = 'Images/'
        else:
            if not imgdir.endswith('/'):
                imgdir = imgdir + '/'

        # Define path to processing directories
        dirs = []
        for day in days:
            procdir = os.path.join(monthdir, day, imgdir)
            # Check full image path
            if not os.path.isdir(procdir):
                print('\nSkipping non-existent directory {}'.format(procdir))
                continue
            else:
                dirs.append(procdir)

    # Make sure there is at least one directory to process
    if len(dirs) < 1:
        raise ConfigError('No valid directories found.')

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

    # Check if smear time is given
    if setup['smear time'] is None:
        setup['smear time'] = 900. #[s]

    # Check list of sky catalogs
    catflag = False
    if stages['catalog matching']:
        # all available catalogs
        catalog_opts = sorted(catalogio.catalog_dict.keys())
        # If catalogs = [], use all of them
        if len(setup['catalogs']) < 1:
            setup['catalogs'] = catalog_opts
            catflag = False
        else:
            catflag = True
            # Make sure requested catalogs exist
            try:
                for cat in setup['catalogs']:
                    if type(cat) != str:
                        cat = str(cat)
                    cat = cat.lower()
                    if cat not in catalog_opts:
                        print('\nCurrently available catalogs: {}\n'.
                              format(catalog_opts))
                        raise ConfigError('Catalog {} is not a valid option'.
                                          format(cat))
            except TypeError:
                raise ConfigError('Please provide a list of valid sky '
                                  'catalogs.')

    # Check required source finding parameters
    if sfparams['mode'] != 'default' and sfparams['mode'] != 'minimize_islands':
        raise ConfigError('Source finding mode must be default or '
                          'minimize_islands.')
    if sfparams['scale'] < 0 or sfparams['scale'] > 10:
        raise ConfigError('The image radius scale factor must be a number '
                          'between 0 and 10.')
    if sfparams['borderpad'] is None or sfparams['borderpad'] < 0:
        sfparams['borderpad'] = 3 #default to 3 pixels

    # Check if beam corrected option is set
    if opts['beam corrected'] is None:
        opts['beam corrected'] = False

    # Check if always associate option is set
    if opts['always associate'] is None:
        opts['always associate'] = False

    # Check if save beam image option is set
    if opts['save beam image'] is None:
        opts['save beam image'] = False

    # Set default QA requirements if not specified
    if opts['quality checks']:
        if qaparams['min nvis'] is None:
            qaparams['min nvis'] = 1300.
        else:
            try:
                qaparams['min nvis'] = float(
                    qaparams['min nvis'])
            except ValueError:
                raise ConfigError('min nvis must be a number.')
        if qaparams['max sensitivity metric'] is None:
            qaparams['max sensitivity metric'] = 3000.
        else:
            try:
                qaparams['max sensitivity metric'] = float(
                    qaparams['max sensitivity metric'])
            except ValueError:
                raise ConfigError('max sensitivity metric must be a number.')
        if qaparams['max beam axis ratio'] is None:
            qaparams['max beam axis ratio'] = 4.
        else:
            try:
                qaparams['max beam axis ratio'] = float(
                    qaparams['max beam axis ratio'])
            except ValueError:
                raise ConfigError('max beam axis ratio must be a number.')
        if qaparams['max source count metric'] is None:
            qaparams['max source count metric'] = 10.
        else:
            try:
                qaparams['max source count metric'] = float(
                    qaparams['max source count metric'])
            except ValueError:
                raise ConfigError('max source count metric must be a number.')
        if qaparams['min niter'] is None:
            qaparams['min niter'] = 1000.
        else:
            try:
                qaparams['min niter'] = float(
                    qaparams['min niter'])
            except ValueError:
                raise ConfigError('min niter must be a number.')
        if qaparams['min bpix'] is None:
            qaparams['min bpix'] = 2.8
        else:
            try:
                qaparams['min bpix'] = float(qaparams['min bpix'])
            except ValueError:
                raise ConfigError('min bpix must be a number.')
        if qaparams['max bpix'] is None:
            qaparams['max bpix'] = 7
        else:
            try:
                qaparams['max bpix'] = float(qaparams['max bpix'])
            except ValueError:
                raise ConfigError('max bpix must be a number.')

    return stages, opts, setup, sfparams, qaparams, dirs, catflag


def dbinit(dbname, user, overwrite, qaparams, safe_override=False):
    """Creates a psycopg2 connection object to communicate
    with the PostgreSQL database. If no database with the
    provided name exists, the user is prompted to create a
    new one. The "radcat" schema which holds all the radio
    catalogs in tables is created at this stage if it does not
    already exist. The user will be prompted to verify deletion
    of all current tables if the database exist and the *overwrite*
    option in the configuration file is ``True``. All necessary
    tables, functions, and triggers are created through a call to 
    ``database.createdb.create()``.

    Parameters
    ----------
    dbname : str
        Name of the PostgreSQL database.
    user : str
        Username for the PostgreSQL database connection.
    overwrite : bool
        If ``True``, tables and data will be deleted and re-created,
        assuming there is a pre-existing database of name "dbname".
        If ``False``, the existing database with "dbname" is used.
    qaparams : dict
        The dictionary of image quality requirements specified by the user
        in the configuration file. The values are written to the
        database **error** table.
    safe_override : bool, optional
        If ``True``, this overrides the 'safe' parameter in
        ``database.createdb.create()``. Default value is ``False``.

    Returns
    -------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    """
    try:
        # DB exists
        conn = psycopg2.connect(host='localhost', database=dbname, user=user)
        logger.info('Connected to database {}.'.format(dbname))
        if not overwrite:
            logger.info('Using existing database {}.'.format(dbname))
        else:
            logger.info('Overwriting existing database tables.')
            if safe_override:
                createdb.create(conn, qaparams, safe=True)
            else:
                # This will prompt warning in create function
                createdb.create(conn, qaparams, safe=False)
        # Check for sky catalogs by verifying schema exists
        cur = conn.cursor()
        cur.execute('''SELECT EXISTS(SELECT 1 FROM pg_namespace 
            WHERE nspname = 'radcat');''')
        if not cur.fetchone()[0]:
            cur.close()
            logger.info('Radio catalog schema "radcat" not found. '
                        'Creating tables now...')
            radcatdb.create(conn)
        else:
            cur.close()
    except psycopg2.OperationalError:
        # DB does not yet exist
        if safe_override:
            makenew = 'yes'
        else:
            makenew = input('\nCreate new database {}? '.format(dbname))
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
            logger.info('Connected to new database {}.'.format(dbname))
            createdb.create(conn, qaparams, safe=True)
            logger.info('Adding radio catalogs to "radcat" schema...')
            radcatdb.create(conn)
        else:
            logger.info('No new database created.')
            raise ConfigError('Cannot access database {}'.format(dbname))

    return conn


def vlite_unique(conn, src, image_id, radius):
    """This function adds a VLITE unique (VU) source, or a source
    with no other radio catalog matches, to the **vlite_unique** table. 
    After adding a new VU source, the **image** table is queried
    to find all previously processed images in which the VU source
    could have been detected based on field-of-view size. The 
    **vlite_unique** table, therefore, keeps a record
    of every image in which a VU source was in the field-of-view
    and whether it was detected.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    src : ``database.dbclasses.DetectedSource`` instance
        The VU source from the **assoc_source** table.
    image_id : int
        Id number of the image in which the VU source
        was detected.
    radius : float
        Radius (in degrees) of the image's field-of-view. Used in
        querying the **image** table.

    Returns
    -------
    src : ``database.dbclasses.DetectedSource`` instance
        The same VU src with updated 'detected' attribute.
    """
    # Start by checking if the VU source is already in the VU table
    existing = dbio.check_vlite_unique(conn, src.id)
    if not existing:  # returned empty list
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
            if not entry[0][2]:  # if not detected
                # Update entry to detected = True
                src.detected = True
                dbio.add_vlite_unique(conn, src, image_id, update=True)

    return src


def iminit(conn, imobj, save, qa, qaparams, reproc, stages, scale, nside, skymap):
    """This function handles ingestion of the image metadata
    into the database **image** table and represents the first
    stage in the VLITE Database Pipeline.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
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
    scale : float
        Fraction between 0 and 1 of the image radius to use.
        The full size of the image field-of-view is multiplied
        by this number.
    nside : int
        Nside parameter of healpy GSM skymap
    skymap : float array
        Healpy format GSM map for setting tsky of image table

    Returns
    -------
    imobj : ``database.dbclasses.Image`` instance
        Image object with `id` attribute updated after insertion into
        the database **image** table, or ``None`` if the image
        has already been processed and will not be reprocessed.
    """
    # STAGE 1 -- Add image to table, 1st quality check
    logger.info('**********************')
    logger.info('STAGE 1: READING IMAGE')
    logger.info('**********************')

    # Set tsky
    imobj.set_tsky(nside, skymap)

    # Set the radius size
    imobj.set_radius(scale)

    # Run image quality checks
    if qa:
        imobj.image_qa(qaparams)
    else:
        imobj.error_id = None
        pass

    # SHOULD THIS HAPPEN FIRST????
    # Is the image in the database?
    status = dbio.status_check(conn, imobj.filename)

    # Only adding image to DB - add or update w/o deleting sources
    if not any(stages.values()):
        if status is None:
            # image not in DB, so add it
            imobj = dbio.add_image(conn, imobj, status)  # branch 2
            global branch
            branch = 2
        else:
            if reproc:
                if save:
                    # already processed, but re-doing -- sources NOT deleted
                    imobj = dbio.add_image(conn, imobj, status)  # branch 4
                else:
                    # just initialize if not writing to DB
                    logger.info('Initializing image.')
                branch = 4
            else:
                # already processed & not re-doing
                logger.info('Image already in database. Moving on...')
                imobj = None  # branch 3
                branch = 3
    else:  # Running at least one stage
        if status is None:
            if stages['source finding']:
                if save:
                    # image not in DB; planning to source find & write to DB
                    imobj = dbio.add_image(conn, imobj, status)
                else:
                    # just initialize if not writing to DB
                    logger.info('Initializing image.')
                branch = 6
            else:
                # image not in DB, but not running source finding --> quit
                logger.error('ERROR: Image {} not yet processed. '
                             'Source finding must be run before other stages.'.
                             format(imobj.filename))
                imobj = None  # branch 5
                branch = 5
        else:
            if stages['source finding']:
                if status[1] == 1:
                    # image in DB, but no SF yet
                    if save:
                        imobj = dbio.add_image(conn, imobj, status)
                    else:
                        logger.info('Initializing image.')
                    branch = 8
                else:
                    if reproc:
                        if save:
                            # image in DB, delete old SF results
                            imobj = dbio.add_image(conn, imobj, status,
                                                   delete=True)
                        else:
                            # just initialize if not writing to DB
                            logger.info('Initializing image.')
                        branch = 8
                    else:
                        # image has SF results & not re-processing
                        logger.info('Image already processed. Moving on...')
                        imobj = None  # branch 9
                        branch = 9
            else:
                # not running SF -- stage must be > 1
                if status[1] > 1:
                    logger.info('Initializing image.')
                    imobj.id = status[0]
                    imobj.stage = status[1]
                    imobj.radius = status[2]
                else:
                    logger.error('ERROR: Image {} does not have sources '
                                 'extracted yet. Source finding must be run '
                                 'before other stages.'.format(imobj.filename))
                    imobj = None  # branch 7
                    branch = 7

    # Stop if the image failed a quality check except 6
    if qa:
        if imobj is not None and imobj.error_id is not None:
            if imobj.error_id != 6:
                imobj = None

    return imobj


def srcfind(conn, imobj, sfparams, save, qa, qaparams, opts, pbdic, beamsdir):
    """Runs PyBDSF source finding and inserts source 
    fit parameters into database  **detected_source**,
    **detected_island**, and **corrected_flux** tables,
    if the *save to database* option is set to ``True``. 
    This function represents the second stage of the
    VLITE Database Pipeline.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info.
    sfparams : dict
        Specifies any non-default PyBDSF parameters to be 
        used in source finding.
    save : bool
        If ``True``, the source fit parameters are written and saved
        to the database **detected_island**, **detected_source**, and
        **corrected_flux** tables. The PyBDSF files are always 
        written out.
    qa : bool
        If ``True``, quality checks are run on the source finding
        results.
    qaparams : dict
        User-specified requirements from the configuration file for
        the source finding quality checks.
    opts : dict
        User-specified options from the configuration file. Needed
        to decide on applying primary beam correction
    pbdic : dict
ass        Dictionary of ``sourcefinding.beam_tools.pribeam`` instances
        with the primary beams
    beamsdir : str
        Directory of .../Beams/ dir for primary beam images
    Returns
    -------
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with updated attributes
        from the source finding results.
    sources : list
        List of ``database.dbclasses.DetectedSource`` objects.
        Attributes of each object are set from the PyBDSF
        output object.
    """
    # STAGE 2 -- Source finding + 2nd quality check
    logger.info('***********************')
    logger.info('STAGE 2: SOURCE FINDING')
    logger.info('***********************')

    # Initialize source finding image object
    bdsfim = runbdsf.BDSFImage(imobj.filename, **sfparams)
    # Run PyBDSF source finding
    if sfparams['mode'] == 'minimize_islands':
        out = bdsfim.minimize_islands()
    else:
        out = bdsfim.find_sources()

    # Update stage
    imobj.stage = 2

    if out is not None:
        # Write PyBDSF files to daily directory
        runbdsf.write_sources(out)
        # Translate PyBDSF output to DetectedSource objects
        sources = dbclasses.translate(imobj, out)
        # Drop sources outside the (scaled) image FOV (radius)
        # Why did PyBDSF give sources with 0 positional errors?
        #  Filter these out. How common is this? Are we losing real source
        #  detections? Was this only older PyBDSF versions?
        if imobj.filename.endswith('IPln1.fits'):
            sources = [src for src in sources if
                       src.dist_from_center <= imobj.radius and
                       src.e_ra > 1e-7 and src.e_dec > 1e-7 and
                       src.xpix > sfparams['borderpad'] and src.ypix > sfparams['borderpad'] and
                       src.xpix < imobj.naxis1 - sfparams['borderpad'] and
                       src.ypix < imobj.naxis2 - sfparams['borderpad']]
            logger.info(' -- {}/{} sources are inside the circular FOV '
                        'with radius {} degree(s)'.format(
                            len(sources), out.nsrc, imobj.radius))
        # Add PyBDSF defined attributes to Image object
        imobj.rms_box = str(out.rms_box)
        imobj.nsrc = len(sources)
        # Run quality checks, part 2
        if qa:
            imobj.source_qa(sources, qaparams)
    else:
        # PyBDSF failed to process
        sources = None
        imobj.error_id = 7

    # Stop if the image failed the source count QA
    if qa:
        if imobj.error_id is not None:
            if imobj.error_id != 6:
                sources = None
    else:
        imobj.error_id = None

    # Determine which sources were CLEANed
    if sources is not None:
        logger.info('Checking which sources were CLEANed.')
        imobj, sources = radioxmatch.check_clean(conn, sources, imobj)

    # Determine nearest neighbors
    if sources is not None:
        logger.info('Calculating nearest neighbors.')
        sources = radioxmatch.nearestneigh(sources)

    if save:
        # Add source fit parameters to database tables
        dbio.add_sources(conn, imobj, sources)
        if sources is not None:
            # Compute beam corrected fluxes, compactness,
            #  & write to corrected_flux table
            # To beam correct or not to beam correct...
            if opts['beam corrected']:  # images already are
                logger.info('Image already beam corrected.')
            else:
                #update pb_flag in image table
                dbio.update_pbflag(conn, imobj)
                #update ass_flag?
                if imobj.pb_flag == False and imobj.ass_flag:
                    if opts['always associate'] == False:
                        imobj.ass_flag = False
                        dbio.update_assflag(conn, imobj)
                if imobj.pb_flag: #if flag True, ok to calc primary beam image
                    logger.info('Calculating beam image and correcting all ' 
                            'flux measurements for primary beam response.')
                    imobj.set_beam_image(pbdic)
                    #create and write fits file 
                    if opts['save beam image']:
                        logger.info('Writing primary beam image file')
                        hdu2 = fits.open(imobj.filename, mode='readonly')
                        hdr2 = hdu2[0].header
                        hdrbm = fits.Header()
                        hdrbm['CTYPE1'] = hdr2['CTYPE1']
                        hdrbm['CDELT1'] = hdr2['CDELT1']
                        hdrbm['CRPIX1'] = hdr2['CRPIX1']
                        hdrbm['CROTA1'] = hdr2['CROTA1']
                        hdrbm['CRVAL1'] = hdr2['CRVAL1']
                        hdrbm['CTYPE2'] = hdr2['CTYPE2']
                        hdrbm['CDELT2'] = hdr2['CDELT2']
                        hdrbm['CRPIX2'] = hdr2['CRPIX2']
                        hdrbm['CROTA2'] = hdr2['CROTA2']
                        hdrbm['CRVAL2'] = hdr2['CRVAL2']
                        hdrbm['INSTRUME'] = 'VLITE'
                        hdrbm['OBSERVER'] = 'VLITE'
                        hdrbm['DATE-OBS'] = imobj.obs_date
                        hdrbm['DATE-MAP'] = imobj.map_date
                        hdrbm['EPOCH'] = hdr2['EPOCH']
                        hdrbm['EQUINOX'] = hdr2['EQUINOX']
                        hdrbm['GLON'] = imobj.glon
                        hdrbm['GLAT'] = imobj.glat
                        hdrbm['CONFIG'] = imobj.config
                        hdrbm['PRIBAND'] = imobj.priband
                        hdrbm['MJDTIME'] = imobj.mjdtime
                        hdrbm['TAU_TIME'] = imobj.tau_time
                        hdrbm['DURATION'] = imobj.duration
                        hdrbm['HISTORY'] = 'Max beam pixel = %e' % np.max(imobj.bmimg)
                        hdrbm['HISTORY'] = 'smear time = %s sec' % str(imobj.smeartime)
                        hdrbm['HISTORY'] = 'nbeam = %d' % imobj.nbeam
                        pbp = '[' + str(imobj.pbparangs[0])
                        pbw = '[' + str(imobj.pbweights[0])
                        pbz = '[' + str(imobj.pbza[0])
                        for ind in range(1,imobj.nbeam):
                            pbp += (',' + str(imobj.pbparangs[ind]))
                            pbw += (',' + str(imobj.pbweights[ind]))
                            pbz += (',' + str(imobj.pbza[ind]))
                        pbp += ']'
                        pbw += ']'
                        pbz += ']'
                        hdrbm['HISTORY'] = 'pbparangs = %s' % pbp
                        hdrbm['HISTORY'] = 'pbweights = %s' % pbw
                        hdrbm['HISTORY'] = 'pbza = %s' % pbz
                        if imobj.ninterval is not None:
                            hdrbm['HISTORY'] = 'ninterval = %d' % imobj.ninterval
                        if imobj.max_dt is not None:
                            hdrbm['HISTORY'] = 'max_dt = %s' % str(imobj.max_dt)
                        hdrbm['HISTORY'] = 'Generated with vdp.py version %s ' % __version__
                        if imobj.vcss:
                            hdrbm['HISTORY'] = 'This is a VCSS snapshot primary beam image.'
                        hdubm = fits.PrimaryHDU(np.float32(imobj.bmimg),header=hdrbm)
                        hdulbm = fits.HDUList([hdubm])
                        btmp = imobj.filename.split('/')[-1].split('.')
                        bnametmp = ''
                        for i in range(len(btmp)-2):
                            bnametmp += (btmp[i]+'.')
                        beamname = beamsdir+bnametmp+'beam.fits'
                        hdulbm.writeto(beamname,overwrite=True)
                        hdu2.close()
                        hdulbm.close()
                    #update pribeam image values in image table
                    dbio.update_pbvalues(conn, imobj)
                    for src in sources:
                        src.correct_flux(imobj,pbdic)
                    # Check if beam solid angle good for association
                    # Don't check VCSS
                    if not imobj.vcss:
                        aveza = np.sum(np.array(imobj.pbza)*np.array(imobj.pbweights))
                        bsa = 1.1331*imobj.bmin*imobj.bmaj*np.cos(np.radians(aveza)) \
                            / dbclasses.res_dict[imobj.config][imobj.cycle]['bsanorm']
                        if bsa < 0.78 or bsa > 1.75:
                            if not opts['always associate']:
                                imobj.ass_flag = False
                                dbio.update_assflag(conn, imobj)
                        else:
                            imobj.ass_flag = True
                            dbio.update_assflag(conn, imobj)
                else:
                    logger.info('Image pb_flag False, cannot calculate beam image.')
            # Calc SNR, compactness, add to table. These don't require beam correction
            for src in sources:
                src.calc_snr()
                src.calc_compactness(imobj)
                dbio.add_corrected(conn, src)

    return imobj, sources


def srcassoc(conn, imobj, sources, save, sfparams):
    """Associates through positional cross-matching sources
    extracted from the current image with previously detected
    VLITE sources stored in the **assoc_source** database table.
    This function represents the third stage of the VLITE
    Database Pipeline.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    sources : list
        List of ``database.dbclasses.DetectedSource`` objects.
        Attributes of each object are from the PyBDSF fit results.
    save : bool
        If ``True``, the **assoc_source** table is updated with the
        association results and the 'assoc_id' is updated in the
        **detected_source** table. If ``False``, no results are
        saved to the database.
    sfparams : dict
        Specifies any non-default PyBDSF parameters to be 
        used in source finding. Passed to nullfind

    Returns
    -------
    detected_unmatched : list
        List of new VLITE detected sources. 
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with updated `stage` attribute.
    """
    # STAGE 3 -- Source association
    logger.info('***************************')
    logger.info('STAGE 3: SOURCE ASSOCIATION')
    logger.info('***************************')

    # If no sources just return
    if len(sources)==0:
        logger.info('No sources to associate!')
        # Update stage in image table
        imobj.stage = 3
        dbio.update_stage(conn, imobj)
        return sources, imobj

    # Limit cone search radius to image FOV or 3 deg for VLITE
    #  (because beam corrections only reliable to 3 deg)
    # Use bigger cone search radius for VCSS mosaics
    if imobj.filename.endswith('IPln1.fits'):
        if imobj.radius < 3:
            radius = imobj.radius  # deg
        else:
            radius = 3. # deg
    else: # mosaics:
        radius = 3.  # deg


    # Associate current sources with existing VLITE catalog
    detected_matched, detected_unmatched, assoc_matched, assoc_unmatched \
        = radioxmatch.associate(conn, sources, imobj, radius, save)

    if save:
        # Update assoc_id col for matched detected sources & corrected_flux
        if detected_matched:
            dbio.update_detected_associd(conn, detected_matched)
            dbio.update_corrected_associd(conn, detected_matched)
        # Add new (unmatched) detected sources to assoc_source table
        if detected_unmatched:
            # Updates assoc_id attribute
            detected_unmatched = dbio.add_assoc(conn, detected_unmatched)
        # Update matched assoc_source positions, etc.
        if assoc_matched:
            dbio.update_matched_assoc(conn, assoc_matched)
            # update associated source light curve metrics
            dbio.update_lcmetrics(conn, assoc_matched)
        # Check for null detections
        # if assoc_unmatched:
        #     nullfind(conn, imobj, sfparams, save, assoc_unmatched)
        # Check for VLITE unique (VU) sources that weren't detected in image
        # **************************************************
        # Nov-2021: Removing this part of the code. Rebuilding VCSS
        #   snapshot assoc_source tables messes up VU table. VU table
        #   not used anyway
        '''
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
        '''
        # **************************************************
        # Update stage in image table
        imobj.stage = 3
        dbio.update_stage(conn, imobj)

    return detected_unmatched, imobj


def catmatch(conn, imobj, sources, catalogs, save, catflag):
    """Performs positional cross-matching of VLITE 
    detected sources to other radio sky survey catalogs.
    This function represents the fourth and final stage
    in the VLITE Database Pipeline.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info.
    sources : list
        VLITE detected sources to be matched to other
        radio catalog sources.
    catalogs : list
        Names of the radio catalogs to use.
    save : bool
        If ``True``, match results are recorded in the
        **catalog_match** table and the **assoc_source** table is
        updated. VLITE unique sources with no sky catalog match are
        inserted into the **vlite_unique** table. If ``False``,
        results are printed to the terminal and no changes are
        made to the database.
    catflag : bool
        If True, match sources to user-given catalogs
        regardless of resolution.
    """
    # STAGE 4 -- Sky catalog cross-matching
    logger.info('*********************************')
    logger.info('STAGE 4: MATCHING TO SKY CATALOGS')
    logger.info('*********************************')

    catalogs = [catalog.lower() for catalog in catalogs]
    if catflag:
        filtered_catalogs = catalogs
    else:
        # Filter catalogs by resolution
        filtered_catalogs = radioxmatch.filter_catalogs(
            conn, catalogs, imobj)
    # Remove catalogs that have already been checked for this image
    if save:
        new_catalogs = dbio.update_checked_catalogs(
            conn, imobj.id, filtered_catalogs)
    else:
        new_catalogs = catalogs
    if not new_catalogs:
        logger.info('All specified catalogs with appropriate resolution '
                    'have already been checked for matches.')
        if save:
            imobj.stage = 4
            dbio.update_stage(conn, imobj)
        return

    logger.info('Using the following catalogs for cross-matching: {}'.format(
        new_catalogs))

    if not sources:
        logger.info('No new VLITE sources to match.')
        if save:
            imobj.stage = 4
            dbio.update_stage(conn, imobj)
        return

    # Limit cone search radius to image FOV for VLITE, bigger for VCSS mosaics
    if imobj.filename.endswith('IPln1.fits'):
        if imobj.radius < 3:
            radius = imobj.radius  # deg
        else:
            radius = 3. # deg
    else: # mosaics:
        radius = 3.  # deg
    #if imobj.filename.endswith('IPln1.fits'):
    #    radius = imobj.radius  # deg
    #else:
    #    radius = 3.  # deg

    # Cross-match VLITE sources with each catalog
    for catalog in new_catalogs:
        try:
            sources, catalog_matched = radioxmatch.catalogmatch(
                conn, sources, catalog, imobj, radius, save)
        except TypeError:
            # No sky catalog sources extracted, move on to next catalog
            continue
        if save:
            # Add results to catalog_match table
            dbio.add_catalog_match(conn, catalog_matched)
    if save:
        # Update assoc_source nmatches
        dbio.update_assoc_nmatches(conn, sources)
        # **************************************************
        # Nov-2021: Removing this part of the code. Rebuilding VCSS
        #   snapshot assoc_source tables messes up VU table. VU table
        #   not used anyway
        '''
        # Check for new VLITE unique (VU) sources from this image
        for src in sources:
            if src.nmatches == 0:
                src = vlite_unique(conn, src, imobj.id, imobj.radius)
        '''
        # **************************************************
        # Update stage
        imobj.stage = 4
        dbio.update_stage(conn, imobj)

    return


def nullfind(conn, imobj, sfparams, save, asrcs):
    """Runs PyBDSF source finding in force-fitting mode
    at locations of unmatched associated sources and 
    inserts null detections into **detected_null** table 
    if the *save to database* option is set to ``True``. 

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info.
    sfparams : dict
        Specifies any non-default PyBDSF parameters to be 
        used in source finding. Will pass a copy with src_ra_dec,
        which gives the force-fit coordinates, to process_image().
    save : bool
        If ``True``, the missed detections are saved
        to the database **detected_null** tables.
    asrcs : list
        List of ``database.dbclasses.DetectedSource`` objects,
        the unmatched associated sources to force-fit.
    """

    logger.info('***********************')
    logger.info('CHECKING %d FOR NULLS' % len(asrcs))
    logger.info('***********************')

    # Add the PyBDSF ``src_ra_dec'' parameter to sfparams as
    # list of tuples with (RA,Dec) coords for force-fitting.
    coords = []
    for src in asrcs:
        coords.append((src.ra, src.dec))

    nfparams = dict(sfparams)  # copy params dict for null finding
    nfparams['src_ra_dec'] = coords
    # nfparams['fix_to_beam'] = True # force Gaussians to beam size
    # set radius of aperture for aperture flux and error calculation
    nfparams['aperture'] = sqrt(imobj.bmin*imobj.bmaj)/imobj.pixel_scale

    # Initialize source finding image object
    bdsfim = runbdsf.BDSFImage(imobj.filename, **nfparams)
    # Run PyBDSF source finding
    out = bdsfim.find_sources()

    nulls = None
    if out is not None:
        nulls = []
        # Translate PyBDSF output to DetectedSource objects
        sources = dbclasses.translate_null(imobj, out, coords)
        # Beam correct fluxes
        logger.info('Correcting null flux measurements for primary beam '
                    'response.')
        for n, src in enumerate(sources):  # *should be* in order with asrcs
            src.correct_flux_null(imobj.pri_freq)
            # Calculate pseudo-SNRs
            if src.total_flux < 0.:
                src.total_flux = 0.
            src.snr = (asrcs[n].ave_total-src.total_flux) / \
                sqrt(asrcs[n].e_ave_total**2 + src.e_total_flux**2)
            # Check for null detections
            if src.snr > 10.0:  # should have been detected
                nulls.append(src)
                # set attributes not set by translate & correct_flux:
                nulls[-1].assoc_id = asrcs[n].id
            # print '%d %f %f  %f %f  %f %f  %.2f  %f %f' % (n,coords[n][0],coords[n][1],src.ra,src.dec,src.total_flux,src.e_total_flux,src.snr,src.dist_from_center,src.polar_angle)

    else:
        # PyBDSF failed to process
        logger.info('PyBDSF failed to process force-fitting!')

    if save:
        if nulls is not None:
            # Add null parameters to database table
            dbio.add_nulls(conn, nulls)

    return


def process(conn, stages, opts, dirs, files, catalogs, sfparams, qaparams, setup, catflag):
    """This function handles the logic and transitions between
    processing stages.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    stages : dict
        Keys are the processing stages (source finding, source assocation,
        and catalog matching) and values are boolean ``True`` or ``False``.
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
        Specifies any non-default PyBDSF parameters to be used in source
        finding.
    qaparams : dict
        User-specified quality requirements or default values defined
        and set in ``cfgparse``.
    setup : dict
        Keys are the setup parameters (root directory, year, month, day,
        files, database name, database user, smear_time and catalogs) and values
        are the user-supplied inputs.
    catflag : boolean
        If True, match sources to all user-given catalogs 
        regardless of resolution.
    """
    global branch
    # Define booleans from stages & opts dictionaries
    sf = stages['source finding']
    sa = stages['source association']
    cm = stages['catalog matching']
    save = opts['save to database']
    qa = opts['quality checks']
    reproc = opts['reprocess']
    rematch = opts['redo match']
    updatematch = opts['update match']
    alwaysass = opts['always associate']

    # Create and enable break handler
    bh = BreakHandler()
    bh.enable()

    # Read GSM map for setting tsky in image table
    fin = open('341.GSM', 'r')
    skymap = np.loadtxt(fin)
    fin.close()
    nside = hp.get_nside(skymap)

    #Create primary beams dictionary
    pbdic = beam_tools.Get_Primary_Beams()

    # Begin loop through daily directories
    i = 0
    for imgdir in dirs:
        # Check if there was a break in image loop
        if bh.trapped:
            logger.info('Pipeline terminated (keyboard interrupt).')
            break
        # Define/make directory for PyBDSF output
        if setup['year'] is None and setup['month'] is None:
            daydir = imgdir
        else:
            daydir = os.path.abspath(os.path.join(imgdir, '..'))
        pybdsfdir = os.path.join(daydir, 'PyBDSF/')
        if not os.path.isdir(pybdsfdir):
            os.system('mkdir '+pybdsfdir)
        # Define/make directory for primary beam images
        beamsdir = os.path.join(daydir, 'Beams/')
        if not os.path.isdir(beamsdir):
            os.system('mkdir '+beamsdir)    

        if not files[0]:
            # Select all images that end with 'IPln1.fits' or 'pbcor.fits' or 'ITime.fits'...
            imglist = [f for f in os.listdir(imgdir) if
                       f.endswith('IPln1.fits') or f.endswith('pbcor.fits') or f.endswith('ITime.fits')]
            # ...or 'IMSC.fits' for the VCSS mosaics
            if len(imglist) < 1:
                imglist = [f for f in os.listdir(imgdir) if
                           f.endswith('IMSC.fits')]
        else:
            imglist = [f for f in files[i]]
        i += 1

        # Make list of objects containing just the image name & mjdtime
        imgmjdlist=[]
        for img in imglist:
            #print(imgdir+img)
            impath = os.path.join(imgdir, img)
            imgmjdlist.append(dbclasses.getimgmjd(impath))

        # Sort list by mjdtime
        imgmjdlist.sort(key=lambda x: x.mjdtime or 0)

        # Begin loop through time-sorted images
        for img in imgmjdlist:
            print(img.filename)
            imobj = dbclasses.init_image(img.filename, alwaysass, setup['smear time'])
            logger.info('_' * (len(imobj.filename) + 10))
            logger.info('Starting {}.'.format(imobj.filename))
            # STAGE 1 -- Add image to database
            imobj = iminit(conn, imobj, save, qa, qaparams, reproc, stages,
                           sfparams['scale'], nside, skymap)
            # Move on to next image if imobj is None
            if imobj is None:
                continue

            # STAGE 2 -- Source finding
            if sf:
                imobj, sources = srcfind(conn, imobj, sfparams, save, qa, qaparams, opts, pbdic, beamsdir)
                # Copy PyBDSF warnings from their log to ours
                with open(imobj.filename+'.pybdsf.log', 'r') as f:
                    lines = f.readlines()
                warnings = [line.strip()
                            for line in lines if 'WARNING' in line]
                if warnings:
                    logger.info('PyBDSF warnings:')
                    for warning in warnings:
                        logger.info(warning)
                # Move PyBDSF output files to PyBDSF directory
                os.system('mv '+imgdir+'/*pybdsf.log '+pybdsfdir+'.')
                if glob.glob(imgdir+'/*pybdsm*'):
                    os.system('mv '+imgdir+'/*pybdsm* '+pybdsfdir+'.')
                if sources is None:
                    # Image failed to process
                    logger.info('sources is none!')
                    continue
                # STAGE 3 -- Source association
                if sa:
                    if imobj.ass_flag:
                        new_sources, imobj = srcassoc(
                            conn, imobj, sources, save, sfparams)
                        # STAGE 4 -- Sky survey catalog cross-matching
                        if cm:  # sf + sa + cm - branch 12, 15
                            # Cross-match new sources only
                            catmatch(conn, imobj, new_sources, catalogs, save, catflag)
                            if glob.glob(imgdir+'*matches.reg'):
                                os.system(
                                    'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                            logger.info('======================================='
                                        '========================================')
                            logger.info('Completed source finding, association, '
                                        'and sky catalog cross-matching on image')
                            logger.info('{}.'.format(imobj.filename))
                            logger.info('======================================='
                                        '========================================')
                            if branch == 6:
                                branch = 12
                            if branch == 8:
                                branch = 15
                        else:  # sf + sa - branch 11, 14
                            if branch == 6:
                                branch = 11
                            if branch == 8:
                                branch = 14
                            logger.info('======================================='
                                        '========================================')
                            logger.info('Completed source finding and association '
                                        'on image')
                            logger.info('{}.'.format(imobj.filename))
                            logger.info('======================================='
                                        '========================================')
                            continue
                    else:
                        logger.info('Image ass_flag = {}, skipping source association'
                                    '/catalog matching'.format(imobj.ass_flag))
                else:
                    if cm:  # sf + cm - branch 10, 13
                        if imobj.ass_flag:
                            catmatch(conn, imobj, sources, catalogs, False, catflag)
                            if glob.glob(imgdir+'*matches.reg'):
                                os.system(
                                    'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                            logger.info('======================================='
                                        '========================================')
                            logger.info('Completed source finding and sky catalog '
                                        'cross-matching to the extracted sources '
                                        'from image')
                            logger.info('{}.'.format(imobj.filename))
                            logger.info('======================================='
                                        '========================================')
                            if branch == 6:
                                branch = 10
                            if branch == 8:
                                branch = 13
                        else:
                            logger.info('Image ass_flag = {}, skipping catalog'
                                        'matching'.format(imobj.ass_flag))
                    else:  # sf only - branch 6, 8
                        logger.info('======================================='
                                    '========================================')
                        logger.info('Completed source finding on image')
                        logger.info('{}.'.format(imobj.filename))
                        logger.info('======================================='
                                    '========================================')
                        continue
            else:  # no sf
                if sa:
                    if imobj.ass_flag:
                        # Get sources from detected_source table
                        sources = dbio.get_image_sources(conn, imobj.id)
                        # Already caught case of no sf but stage < 2 in iminit
                        if imobj.stage == 2:  # no sa has been run yet
                            new_sources, imobj = srcassoc(conn, imobj,
                                                          sources, save)
                            if cm:  # sa + cm - branch 20
                                # Cross-match new sources only
                                catmatch(conn, imobj, new_sources,
                                         catalogs, save, catflag)
                                if glob.glob(imgdir+'*matches.reg'):
                                    os.system(
                                        'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                                logger.info('====================================='
                                            '====================================='
                                            '=====')
                                logger.info('Completed source association and sky '
                                            'catalog cross-matching to the newly '
                                            'detected sources from image')
                                logger.info('{}.'.format(imobj.filename))
                                logger.info('====================================='
                                            '====================================='
                                            '=====')
                                branch = 20
                            else:  # sa only - branch 19
                                logger.info('====================================='
                                            '====================================='
                                            '=====')
                                logger.info('Completed source association for '
                                            'image')
                                logger.info('{}.'.format(imobj.filename))
                                logger.info('====================================='
                                            '====================================='
                                            '=====')
                                branch = 19
                        else:  # stage > 2
                            logger.info("\nNOTE: {}'s".format(imobj.filename))
                            logger.info('sources have already been associated '
                                        'with the existing VLITE catalog.')
                            if cm:  # cm only - branch 21
                                assoc_sources = dbio.get_associated(
                                    conn, sources)
                                if rematch:
                                    # Delete & redo matching
                                    assoc_sources = dbio.delete_matches(
                                        conn, assoc_sources, imobj.id)
                                else:
                                    if not updatematch:
                                        # Cross-match new/un-matched sources only
                                        assoc_sources = [src for src in
                                                         assoc_sources if
                                                         src.nmatches is None
                                                         or src.nmatches == 0]
                                    else:
                                        # Use all sources if updating
                                        pass
                                catmatch(conn, imobj, assoc_sources,
                                         catalogs, save, catflag)
                                if glob.glob(imgdir+'*matches.reg'):
                                    os.system(
                                        'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                                logger.info('====================================='
                                            '====================================='
                                            '=====')
                                logger.info('Completed sky catalog cross-matching '
                                            'for image')
                                logger.info('{}.'.format(imobj.filename))
                                logger.info('====================================='
                                            '====================================='
                                            '=====')
                                branch = 21
                            else:
                                # branch 18
                                branch = 18
                                continue
                    else:
                        logger.info('Image ass_flag = {}, skipping '
                                    'catalog matching'.format(imobj.ass_flag))
                else:
                    if cm:  # cm only - branch 17
                        pass
                    else:  # branches 2, 4 continued
                        continue
                    if imobj.stage > 2:
                        # Get detected, then assoc sources
                        sources = dbio.get_image_sources(conn, imobj.id)
                        assoc_sources = dbio.get_associated(conn, sources)
                        if rematch:
                            # Delete & redo matching
                            assoc_sources = dbio.delete_matches(
                                conn, assoc_sources, imobj.id)
                            branch = 17.2
                        else:
                            if not updatematch:
                                # Cross-match new/un-matched sources only
                                assoc_sources = [src for src in assoc_sources
                                                 if src.nmatches is None or
                                                 src.nmatches == 0]
                                branch = 17.1
                            else:
                                # Use all sources if updating
                                branch = 17.3
                                pass
                        catmatch(conn, imobj, assoc_sources, catalogs, save, catflag)
                        if glob.glob(imgdir+'*matches.reg'):
                            os.system(
                                'mv '+imgdir+'*matches.reg '+pybdsfdir+'.')
                        logger.info('======================================='
                                    '========================================')
                        logger.info('Completed sky catalog cross-matching for '
                                    'image')
                        logger.info('{}.'.format(imobj.filename))
                        logger.info('======================================='
                                    '========================================')
                    else:  # branch 16
                        logger.info('======================================='
                                    '========================================')
                        logger.error('ERROR: Source association must be run '
                                     'before catalog cross-matching for image')
                        logger.error('{}.'.format(imobj.filename))
                        logger.info('======================================='
                                    '========================================')
                        branch = 16
                        continue

            # Check whether there was a break
            if bh.trapped:
                break

    # Disable the break handler
    bh.disable()

    return


def main():
    """One function to rule them all."""
    # Set required & optional command line arguments
    parser = argparse.ArgumentParser(
        description='Run the VLITE Database Pipline (vdp)')
    parser.add_argument('config_file', help='the YAML configuration file')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='stops printing of messages to the console')
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
    parser.add_argument('--unassoc_image', action='store_true',
                        help='unassociates sources in the specified image(s)'
                        'and updates the database')
    parser.add_argument('--manually_add_match', action='store_true',
                        help='manually add catalog matching results for '
                        'VLITE source(s) after follow-up')
    parser.add_argument('--add_catalog', action='store_true',
                        help='adds any new sky survey catalogs to a table in '
                        'the database "radcat" schema')
    #    parser.add_argument('--update_pbcor', action='store_true',
    #                        help='update corrected_flux table with primary '
    #                        'beam corrections')
    args = parser.parse_args()

    # Start the timer
    start_time = datetime.now()

    # Parse run configuration file
    stages, opts, setup, sfparams, qaparams, dirs, catflag = cfgparse(args.config_file)

    # Initialize logger handlers for console & file
    logfile = str(setup['year']) + str(setup['month']).zfill(2) + '.log'
    logpath = os.path.join(setup['root directory'], logfile)
    loggerinit(logpath, args.quiet)
    logger.info('')
    logger.info('#' * (len(logpath) + 10))
    logger.info('Starting the VLITE Database Pipeline.')
    logger.info('Log file: {}'.format(logpath))
    logger.info('#' * (len(logpath) + 10))

    # Find existing/create/overwrite database
    #if any([args.remove_catalog_matches, args.remove_source,
    #        args.remove_image, args.manually_add_match, args.add_catalog, args.update_pbcor]):
    if any([args.remove_catalog_matches, args.remove_source,
            args.remove_image, args.unassoc_image, args.manually_add_match, args.add_catalog]):
        opts['overwrite'] = False
    conn = dbinit(setup['database name'], setup['database user'],
                  opts['overwrite'], qaparams, safe_override=args.ignore_prompt)

    # Option to remove matching results for sky catalogs
    if args.remove_catalog_matches:
        catalogs = input('\nFor which catalogs would you like to remove '
                             'matching results? (List catalogs separated by '
                             'a comma.)\n')
        cat_list = [cat.lower() for cat in catalogs.split(', ')]
        logger.info('Removing matching results for {}...'.format(catalogs))
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
        inp = input('\nPlease enter the id number(s) (i.e. 1, 2, 3) '
                        'of the source(s) you wish to remove from the '
                        'database assoc_source table, or provide a text '
                        'file with one id number per line:\n')
        try:
            asid_list = [int(asid) for asid in inp.strip('[]').split(',')]
        except ValueError:
            with open(inp, 'r') as f:
                text = f.read()
            asid_list = [int(asid) for asid in text.strip().split('\n')]
        logger.info('Removing row(s) {} from the assoc_source table...'.
                    format(asid_list))
        dbio.remove_sources(conn, tuple(asid_list))
        conn.close()
        sys.exit(0)

    # Option to remove images
    if args.remove_image:
        inp = input('\nPlease enter the image(s) filename(s) starting '
                        'at least with the year-month directory (i.e. '
                        '2018-01/15/Images/10GHz.Mrk110.IPln1.fits), or '
                        'provide a text file with one filename per line:\n')
        try:
            images = [re.findall('([0-9]{4}-\S+)', img)[0] for img
                      in inp.split(',')]
        except IndexError:
            print('trying to read ',inp)
            with open(inp, 'r') as f:
                text = f.read()
            images = [img for img in text.strip().split('\n')]
        logger.info('Preparing to remove image(s) {} from the database.'.
                    format(images))
        confirm = input('\nAre you sure? ')
        if confirm == 'y' or confirm == 'yes':
            logger.info('Deleting image(s) from the database...')
            dbio.remove_images(conn, images)
        else:
            logger.info('Doing nothing...')
        conn.close()
        sys.exit(0)

    # Option to unassociate images
    if args.unassoc_image:
        inp = input('\nPlease enter the image(s) filename(s) starting '
                        'at least with the year-month directory (i.e. '
                        '2018-01/15/Images/10GHz.Mrk110.IPln1.fits), or '
                        'provide a text file with one filename per line:\n')
        try:
            images = [re.findall('([0-9]{4}-\S+)', img)[0] for img
                      in inp.split(',')]
        except IndexError:
            print('trying to read ',inp)
            with open(inp, 'r') as f:
                text = f.read()
            images = [img for img in text.strip().split('\n')]
        logger.info('Preparing to unassociate image(s) {} from the database.'.
                    format(images))
        confirm = input('\nAre you sure? ')
        if confirm == 'y' or confirm == 'yes':
            logger.info('Unassociating image(s) from the database...')
            dbio.unassoc_images(conn, images)
        else:
            logger.info('Doing nothing...')
        conn.close()
        sys.exit(0)

    # Option to manually add catalog matching results
    if args.manually_add_match:
        inp = input('\nPlease enter the source assoc_source id, the '
                        'name of the catalog, and, optionally, the id of the '
                        'matched catalog source and the angular separation in '
                        'arcseconds, in that order one per line. NO COMMAS '
                        'Hit "q" when you are done. You may alternatively '
                        'provide a similarly formatted text file with one '
                        'catalog match per line:\n')
        cmatches = []
        while inp != 'q':
            try:
                int(inp.split()[0])
                cmatches.append(inp)
                inp = input()
            except IndexError:
                inp = input()
            except ValueError:
                with open(inp, 'r') as f:
                    text = f.readlines()
                    cmatches = [line.strip('\n') for line in text]
                    break
            if inp == 'q':
                break
        assoc_ids = []
        catalogs = []
        catsrc_ids = []
        separations = []
        for cmatch in cmatches:
            cm = cmatch.split(' ')
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
        logger.info('Adding new catalog matching results for assoc_ids {}...'.
                    format(assoc_ids))
        dbio.update_assoc_nmatches(conn, assoc_ids)
        dbio.add_catalog_match(conn, cmrows)
        conn.close()
        sys.exit(0)

    # Option to add a new sky survey catalog to the database "radcat" schema
    if args.add_catalog:
        radcatdb.create(conn)
        conn.close()
        sys.exit(0)
    '''
    #Needs updating to work with v3.0 beam corrections
    # Option to update corrected_flux table with primary beam corrections
    if args.update_pbcor:
        logger.info('Updating corrected_flux table...')
        dbio.update_corrected(conn,setup)
        conn.close()
        logger.info('Done.')
        sys.exit(0)
    '''
    # Process images
    process(conn, stages, opts, dirs, setup['files'],
            setup['catalogs'], sfparams, qaparams, setup, catflag)

    # Update run_config table & close database connection
    nimages, exec_time = print_run_stats(start_time)
    # Record run configuration parameters
    dbio.record_config(conn, args.config_file, logpath, start_time, exec_time,
                       nimages, stages, opts, setup, sfparams, qaparams)

    conn.close()


if __name__ == '__main__':
    main()
