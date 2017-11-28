"""This is the main script for the VLITE Post-Processing 
Pipeline (P3). It is responsible for reading in the
configuration file, connecting to the the `PostgreSQL` database,
and calling the appropriate processing stages in the correct order.

"""
import os
import sys
import glob
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
    print('\n======================================\n')
    print('Run statistics:')
    # Count the number of images processed
    dbclasses.Image.image_count()
    # Print the runtime
    print('\nTotal runtime: {}\n'.format(datetime.now() - start_time))


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
    opts : tuple of 5 booleans
        Turns on and off options listed in the following order: 
        *save to database, quality checks, overwrite, reprocess, redo match*.
        Items can be ``True`` or ``False``.
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
    res_tol : float
        Resolution tolerance in arcsec. Used to constrain source association
        and catalog matching between sources extracted from images of
        similar spatial resolution.
    """    
    with open(cfgfile, 'r') as stream:
        data = load(stream, Loader=Loader)
        sf = (data['stages'])['source finding']
        sa = (data['stages'])['source association']
        cm = (data['stages'])['catalog matching']
        save = (data['options'])['save to database']
        qa = (data['options'])['quality checks']
        owrite = (data['options'])['overwrite']
        reproc = (data['options'])['reprocess']
        rematch = (data['options'])['redo match']
        rootdir = (data['setup'])['rootdir']
        yr = (data['setup'])['year']
        mo = (data['setup'])['month']
        days = (data['setup'])['day']
        dbname = (data['setup'])['dbname']
        dbusr = (data['setup'])['dbuser']
        catalogs = (data['setup'])['catalogs']
        params = data['pybdsf_params']
        res_tol = (data['matching'])['resolution tolerance']

    stages = (sf, sa, cm)
    opts = (save, qa, owrite, reproc, rematch)

    if not any(stages) and not save:
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
    for stage in stages:
        if isinstance(stage, bool):
            pass
        else:
            raise ConfigError('stage inputs must be True/False or yes/no.')
    for opt in opts:
        if isinstance(opt, bool):
            pass
        else:
            raise ConfigError('option inputs must be True/False or yes/no.')

    # Catch case when no DB is given
    if dbname is None or dbusr is None:
        raise ConfigError('Please provide a database/user name.')

    if cm:
        # Make sure requested catalogs exist
        catalog_opts = catalogio.catalog_list
        try:
            if len(catalogs) < 1:
                catalogs = None
            for cat in catalogs:
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

    # Check res_tol input if running source association or catalog matching
    if sa or cm:
        try:
            # WARNING: yes & no will evaluate to 1.0 & 0.0
            res_tol = float(res_tol)
        except ValueError:
            raise ConfigError('Matching: resolution tolerance must be a '
                              'number.')
        except TypeError: # None -- will assign default value later
            pass
      
    return stages, opts, dirs, dbname, dbusr, catalogs, params, res_tol


def dbinit(dbname, user, overwrite, safe_override=False):
    """Creates a `psycopg2` connection object to communicate
    with the `PostgreSQL` database. If no database with the
    provided name exists, the user is prompted to create a
    new one. This is done by first connecting to the *postgres* 
    database and then calling the SQL CREATE DATABASE command. 
    If the database is new or if overwriting the old one, all 
    necessary tables and triggers are additionally created through 
    a call to `database.createdb.create()`.
    
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
    safe_override : bool
        If ``True``, this overrides the safe boolean. Implemented
        mainly for testing purposes. Default value is ``False``.

    Returns
    -------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    """
    try:
        # DB exists
        conn = psycopg2.connect(host='localhost', database=dbname, user=user)
        print('\nConnected to database {}.'.format(dbname))
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
        if not overwrite:
            print('\nUsing existing database {}.'.format(dbname))
        else:
            print('\nOverwriting existing database tables.')
            if safe_override:
                createdb.create(conn, safe=True)
            else:
                # This will prompt warning in create func
                createdb.create(conn, safe=False)
    except psycopg2.OperationalError:
        # DB does not yet exist
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
            createdb.create(conn, safe=True)
            print('\nCreating new sky catalog tables...')
            skycatdb.create(conn)
        else:
            print('\nNo new database created.\n')
            raise ConfigError('Cannot access database {}'.format(dbname))

    return conn


def iminit(conn, impath, save, qa, reproc, stages):
    """Initializes `database.dbclasses.Image()` object
    and adds row to database image table if the is not
    already in the table or if the image is being reprocessed.
    This function represents the first stage in the 
    Post-Processing Pipeline.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL database connection object.
    impath : str
        Directory path to the fits image file.
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
        Initialized `Image` object with attribute values
        set from header info, or ``None`` if the image
        has already been processed and will not be reprocessed.
    """
    # STAGE 1 -- Initialize image object, add to table, 1st quality check

    # Is the image in the database?
    status = dbio.status_check(conn, impath)

    # Only adding image to DB - add or update w/o deleting sources
    if not any(stages):
        if status is None:
            # image not in DB
            imobj = dbio.add_image(conn, impath, status)
        else:
            if reproc:
                # already processed, but re-doing -- sources NOT deleted
                imobj = dbio.add_image(conn, impath, status)
            else:
                # already processed & not re-doing
                print('\nImage {} already in database. Moving on...'.format(
                    impath))
                imobj = None
    else: # Running at least one stage
        if status is None:
            if stages[0]:
                if save:
                    # image not in DB; planning to source find & write to DB
                    imobj = dbio.add_image(conn, impath, status)
                else:
                    # just initialize if not writing to DB
                    imobj = dbio.init_image(impath)
            else:
                # image not in DB, but not running source finding --> quit
                print('\nERROR: Image {} not yet processed. Source finding '
                      'must be run before other stages.'.format(impath))
                imobj = None
        else:
            if stages[0]:
                if reproc:
                    if save:
                        # image in DB; re-doing SF, so old sources are removed
                        imobj = dbio.add_image(conn, impath,
                                                status, delete=True)
                    else:
                        # just initialize if not writing to DB
                        imobj = dbio.init_image(impath)
                else:
                    # image already in DB & not re-processing
                    print('\nImage {} already processed. Moving on...'.format(
                        impath))
                    imobj = None
            else:
                # not running SF -- stage must be > 1
                if status[1] > 1:
                    imobj = dbio.init_image(impath)
                else:
                    print('\nERROR: Image {} does not have sources extracted '
                          'yet. Source finding must be run before other '
                          'stages.'.format(impath))
                    imobj = None

    # Run quality checks?
    if qa:
        # quality check
        pass
    else:
        pass
    
    return imobj


def srcfind(conn, imobj, params, save, qa):
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
        Initialized `Image` object with attribute values
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
        Initialized `Image` object with attribute values
        set from header info & updated with source finding results.
    sources : list
        List of `database.dbclasses.DetectedSource` objects.
        Attributes of each object are set from the `PyBDSF`
        output object.
    """
    # STAGE 2 -- Source finding + 2nd quality check
    print('\nExtracting sources from {}'.format(imobj.filename))
    # Initialize source finding image object
    bdsfim = runbdsf.BDSFImage(imobj.filename, **params)
    # Run PyBDSF source finding
    if params['mode'] == 'minimize_islands':
        out = bdsfim.minimize_islands()
    else:
        out = bdsfim.find_sources()

    if out is not None:
        # Write PyBDSF(M) files to daily directory
        runbdsf.write_sources(out)
        imobj, sources = dbclasses.translate(imobj, out)
        if save:
            dbio.add_sources(conn, imobj, sources)
        # run more quality checks
        if qa:
            pass
        else:
            pass
    else:
        sources = None

    return imobj, sources


def srcassoc(conn, imobj, sources, res_tol, save):
    """Associates sources extracted from current image with
    previously detected VLITE sources stored in the
    assoc_source database table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Initialized `Image` object with attribute values
        set from header info & updated with source finding results.
    sources : list
        List of `database.dbclasses.DetectedSource` objects.
        Attributes of each object are set from the `PyBDSF`
        output object.
    res_tol : float
        Defines the acceptable range of spatial resolutions in arcsec
        when considering a source from the assoc_source table for
        cross-matching.
    save : bool
        If ``True``, the extracted sources are written and saved
        to the database detected_island and detected_source tables. The
        `PyBDSF` files are always written out.
    
    Returns
    -------
    detected_unmatched : list of `DetectedSource` objects
        Sources extracted from the image which could not be successfully
        associated with previously detected sources. These are added
        to the assoc_table as new detections and then cross-matched
        with other sky survey catalogs. If no match is found, then
        the source becomes a transient candidate.
    imobj : database.dbclasses.Image instance
        Initialized `Image` object with updated stage attribute.
    """
    # STAGE 3 -- Source association
    #search_box = imobj.box_pix2world(imobj.trim_box)
    search_radius = radioxmatch.search_radius(imobj)
    detected_matched, detected_unmatched, assoc_matched, assoc_unmatched \
        = radioxmatch.associate(conn, sources, imobj, search_radius, res_tol)
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
        # Update stage in image table
        imobj.stage = 3
        dbio.update_stage(conn, imobj)

    return detected_unmatched, imobj


def catmatch(conn, imobj, sources, catalogs, res_tol, save):
    # STAGE 4 -- Sky catalog cross-matching
    search_radius = radioxmatch.search_radius(imobj)
    for catalog in catalogs:
        catalog = catalog.lower()
        try:
            sources, catalog_matched = radioxmatch.catalogmatch(
                conn, sources, catalog, imobj, search_radius)
        except TypeError: # no sky catalog sources extracted
            continue
        if save:
            # Add results to catalog_match table
            dbio.add_catalog_match(conn, catalog_matched)
    if save:
        # Update assoc_source.nmatches
        dbio.update_assoc_nmatches(conn, sources)
        # Check for new VLITE unique (VU) sources from this image
        for src in sources:
            if src.nmatches == 0:
                src.detected = True
                # Add source to vlite_unique table
                dbio.add_vlite_unique(conn, src, src.image_id)
                # Find previous images where VU source is in FOV
                src.detected = False
                prev_images = radioxmatch.check_previous(
                    conn, src, search_radius, res_tol)
                for imgid in prev_images:
                    # Ignore current image
                    if imgid[0] != imobj.id:
                        dbio.add_vlite_unique(conn, src, imgid[0])
        # Update stage
        imobj.stage = 4
        dbio.update_stage(conn, imobj)

    return

            
def process(conn, stages, opts, dirs, catalogs, params, res_tol):
    """
    This function handles the data flow logic and transitions
    between processing stages. All individual processing functions
    are called from here. Input parameters come from output of `cfgparse`.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    stages : tuple of 3 booleans
        Indicates which stages to run listed in the following order: 
        *source finding, source association, catalog matching*.
        Items can be ``True`` or ``False``.
    opts : tuple of 5 booleans
        Turns on and off options listed in the following order: 
        *save to database, quality checks, overwrite, reprocess, redo match*.
        Items can be ``True`` or ``False``.
    dirs : list of str
        List of strings specifying paths to daily image directories
        to be processed during the run.
    skycat : str
        Name of the `PostgreSQL` database containing the sky 
        survey catalogs.
    catalogs : list of str
        List of sky survey tables to use when running catalog matching.
    params : dict
        Specifies any non-default `PyBDSF` parameters to be used in source
        finding. 
    res_tol : float
        Resolution tolerance in arcsec. Used to constrain source association
        and catalog matching between sources extracted from images of
        similar spatial resolution.    
    """
    sf, sa, cm = stages
    save, qa, owrite, reproc, rematch = opts
    
    # Begin loop through daily directories
    for imgdir in dirs:
        # Define/make directory for PyBDSF output
        daydir = os.path.abspath(os.path.join(imgdir, '..'))
        pybdsfdir = os.path.join(daydir, 'PyBDSF/')
        if not os.path.isdir(pybdsfdir):
            os.system('mkdir '+pybdsfdir)

        # Select only the images that end with 'IPln1.fits'
        imglist = [f for f in os.listdir(imgdir) if \
                   f.endswith('.0217+738.IPln1.fits')]
        #imglist = [f for f in os.listdir(imgdir) if f.endswith('IPln1.fits')]
        imglist.sort()
        #imglist = imglist[:15]

        # Begin loop through images
        for img in imglist:
            impath = os.path.join(imgdir, img)

            # STAGE 1 -- Initialize image object
            imobj = iminit(conn, impath, save, qa, reproc, stages)
            # Move on to next image if imobj is None
            if imobj is None:
                continue
            
            # STAGE 2 -- Source finding
            if sf:
                imobj, sources = srcfind(conn, imobj, params, save, qa)
                #os.system('rm '+imgdir+'*.crop.fits')
                if sources is None:
                    # Image failed to process
                    with open(pybdsfdir+'failed.txt', 'a') as f:
                        f.write(img+'\n')
                    if save:
                        dbio.pybdsf_fail(conn, imobj)
                    continue
                else:
                    os.system('mv '+imgdir+'*pybdsf.log '+pybdsfdir+'.')
                    if glob.glob(imgdir+'*pybdsm.srl'):
                        os.system('mv '+imgdir+'*pybdsm* '+pybdsfdir+'.')
                # STAGE 3 -- Source association
                if sa:
                    new_sources, imobj = srcassoc(conn, imobj, sources,
                                                  res_tol, save)
                    # STAGE 4 -- Sky survey catalog cross-matching
                    if cm: # sf + sa + cm
                        # Cross-match new sources only
                        catmatch(conn, imobj, new_sources, catalogs,
                                 res_tol, save)
                    else: # sf + sa
                        print('\nSF+SA')
                        continue
                else:
                    # STAGE 4 -- Sky survey catalog cross-matching
                    if cm: # sf + cm
                        if rematch:
                            print('\nSF+CM with rematch')
                            pass
                            #catalogmatch()
                        else:
                            print('\nSF+CM no rematch')
                            pass
                            #catalogmatch() nulls
                    else: # sf only
                        print('\nSF only')
                        continue
            # STAGE 3 -- Source association
            else:
                # no sf, so get sources out of DB
                if sa:
                    # already caught case of no sf but stage < 2
                    if imobj.stage == 2:
                        # run sa for first time on sources
                        pass
                        #sources = some_function()
                        #srcassoc(conn, imobj, sources, res_tol)
                    else: # stage > 2
                        if reproc:
                            # redo sa on sources
                            pass
                            #sources = some_function()
                            #stage3(conn, imobj, sources, res_tol)
                        else:
                            print('\nImage {} already processed. '
                                  'Moving on...'.format(impath))
                            continue
                    # STAGE 4 -- Sky survey catalog cross-matching
                    if cm: # sa + cm
                        #asources = some_other_function()
                        if rematch:
                            print('\nSA+CM with rematch')
                            pass
                            #catalogmatch()
                        else:
                            print('\nSA+CM no rematch')
                            pass
                            #catalogmatch() nulls
                    else: # sa only
                        print('\nSA only')
                        continue
                # STAGE 4 -- Sky survey catalog cross-matching
                else: # cm only
                    if imobj.stage > 2:
                        pass
                        #sources = some_function()
                        #asources = some_other_function()
                        if rematch:
                            print('\nCM only with rematch')
                            pass
                            #catalogmatch()
                        else:
                            print('\nCM only no rematch')
                            pass
                            #catalogmatch() nulls
                    else:
                        print('\nERROR: Source association for image {} must '
                              'be run before catalog cross-matching. '
                              'Alternatively, run source finding first and '
                              'cross-match extracted sources with sky survey '
                              'catalogs.'.format(impath))
                        continue
                        

                '''
                # 4.) Catalog cross-matching
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
    return

            
def main():
    # Start the timer
    start_time = datetime.now()

    try:
        cf = sys.argv[1]
    except IndexError:
        raise ConfigError('Please provide a configuration file.')
    # Parse config file
    stages, opts, dirs, dbname, dbusr, catalogs, params, rtol = cfgparse(cf)

    # Find existing/create/overwrite database
    conn = dbinit(dbname, dbusr, opts[2])

    # Process images
    process(conn, stages, opts, dirs, catalogs, params, rtol)

    conn.close()

    print_run_stats(start_time)


if __name__ == '__main__':
    main()
