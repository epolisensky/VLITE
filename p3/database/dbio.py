"""This module contains methods to prepare and insert
data into the `PostgreSQL` database tables.

"""
import psycopg2
import psycopg2.extras
import dbclasses


def init_image(impath):
    """Initializes an object of the Image class and sets values 
    for its attributes from the fits file header using
    the `header_attrs` class method.

    """
    img = dbclasses.Image(impath)
    data, header = img.read()
    img.header_attrs(header) # Use header info to set attributes

    return img


def status_check(conn, impath):
    """Returns the id, highest completed stage, and radius used
    for source finding from the database image table or ``None``
    if the image filename is not in the database.

    """
    cur = conn.cursor()
    cur.execute('SELECT id, stage, radius FROM Image WHERE filename = %s',
                (impath, ))
    status = cur.fetchone()
    conn.commit()
    cur.close()

    return status


def add_image(conn, impath, status, delete=False):
    """Insert or update rows in database image table.
    
    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    impath : str
        Directory path to the FITS image file.
    status : tuple
        Contains the id and stage for the image's row
        entry in the database image table if it exists.
        Otherwise, status is ``None``.
    delete : bool, optional
        If ``True``, rows in the detected_island database table
        with the appropriate image_id will be deleted,
        cascading to the detected_source table and triggering 
        updates on the assoc_source and catalog_match
        tables. Rows with the same image_id are also deleted 
        from the vlite_unique table. Default value is ``False``.

    Returns
    -------
    img : database.dbclasses.Image instance
        Initialized Image object which was used to
        insert values into the database image table.    
    """
    # Initialize Image object
    img = init_image(impath)    
    
    cur = conn.cursor()
    
    # Add new image to DB
    if status is None:
        print('\nAdding {} to database'.format(img.filename))
        sql = '''INSERT INTO image (
            filename, imsize, obs_ra, obs_dec, pixel_scale, object, obs_date, 
            map_date, obs_freq, primary_freq, bmaj, bmin, bpa, noise, peak, 
            config, nvis, mjdtime, tau_time, duration, radius, nsrc, rms_box, 
            error_id, stage) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id'''
        vals = (img.filename, img.imsize, img.obs_ra, img.obs_dec,
                img.pixel_scale, img.obj, img.obs_date, img.map_date,
                img.obs_freq, img.pri_freq, img.bmaj, img.bmin, img.bpa,
                img.noise, img.peak, img.config, img.nvis, img.mjdtime,
                img.tau_time, img.duration, img.radius, img.nsrc, img.rms_box,
                img.error_id, img.stage)
        cur.execute(sql, vals)
        img.id = cur.fetchone()[0]
    # Update existing image entry
    else:
        print('\nUpdating existing entries for {}'.format(img.filename))
        img.id = status[0]
        sql = '''UPDATE image SET filename = %s, imsize = %s, obs_ra = %s,
            obs_dec = %s, pixel_scale = %s, object = %s, obs_date = %s, 
            map_date = %s, obs_freq = %s, primary_freq = %s, bmaj = %s, 
            bmin = %s, bpa = %s, noise = %s, peak = %s, config = %s, 
            nvis = %s, mjdtime = %s, tau_time = %s, duration = %s, radius = %s,
            nsrc = %s, rms_box = %s, error_id = %s, stage = %s
            WHERE id = %s'''
        vals = (img.filename, img.imsize, img.obs_ra, img.obs_dec,
                img.pixel_scale, img.obj, img.obs_date, img.map_date,
                img.obs_freq, img.pri_freq, img.bmaj, img.bmin, img.bpa,
                img.noise, img.peak, img.config, img.nvis, img.mjdtime,
                img.tau_time, img.duration, img.radius, img.nsrc, img.rms_box,
                img.error_id, img.stage, img.id)
        cur.execute(sql, vals)
        if delete:
            # Delete corresponding sources
            print('\nRemoving previous sources...')
            cur.execute('DELETE FROM detected_island WHERE image_id = %s', (
                img.id, ))

    conn.commit()
    cur.close()

    return img


def add_sources(conn, img, sources):
    """Inserts DetectedSource objects into database
    detected_island and detected_source tables. The Image object
    and image table are also updated with some results
    from the source finding.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.    
    img : database.dbclasses.Image instance
        Initialized Image object which was used to
        insert values into the database image table.
    sources : list of DetectedSource objects
        DetectedSource objects to be inserted into
        the detected_island and detected_source 
        database tables.
    """
    cur = conn.cursor()

    # Update image table -- overwrites all possible updated columns
    sql = '''UPDATE image SET imsize = %s, obs_freq = %s, bmaj = %s, 
        bmin = %s, bpa = %s, noise = %s, radius = %s, nsrc = %s, 
        rms_box = %s, error_id = %s, stage = %s WHERE id = %s'''
    vals = (img.imsize, img.obs_freq, img.bmaj, img.bmin, img.bpa, img.noise,
            img.radius, img.nsrc, img.rms_box, img.error_id, img.stage, img.id)
    cur.execute(sql, vals)
    
    # Insert DetectedSources into detected_source and detected_island tables
    print('\nAdding detected sources to database.')
    for src in sources:
        src.image_id = img.id
        # Insert values into detected_island table
        sql = '''INSERT INTO detected_island (
            isl_id, image_id, total_flux, e_total_flux, 
            rms, mean, resid_rms, resid_mean) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (isl_id, image_id)
            DO NOTHING'''
        vals = (src.isl_id, src.image_id, src.total_flux_isl,
                src.total_flux_islE, src.rms_isl, src.mean_isl,
                src.resid_rms, src.resid_mean)
        cur.execute(sql, vals)

        # Insert values into detected_source table
        sql = '''INSERT INTO detected_source (
            src_id, isl_id, image_id, ra, e_ra, dec, e_dec,
            total_flux, e_total_flux, peak_flux, e_peak_flux, 
            ra_max, e_ra_max, dec_max, e_dec_max, maj, e_maj, 
            min, e_min, pa, e_pa, dc_maj, e_dc_maj, dc_min, e_dc_min,
            dc_pa, e_dc_pa, code, assoc_id) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'''
        vals = (src.src_id, src.isl_id, src.image_id, src.ra, src.e_ra,
                src.dec, src.e_dec, src.total_flux,
                src.e_total_flux, src.peak_flux, src.e_peak_flux,
                src.ra_max, src.e_ra_max, src.dec_max,
                src.e_dec_max, src.maj, src.e_maj, src.min,
                src.e_min, src.pa, src.e_pa, src.dc_maj,
                src.e_dc_maj, src.dc_min, src.e_dc_min, src.dc_pa,
                src.e_dc_pa, src.code, src.assoc_id)
        cur.execute(sql, vals)

    conn.commit()
    cur.close()


def add_assoc(conn, sources):
    """Adds a newly detected VLITE source to the
    assoc_source table and updates the assoc_id
    for that source in the detected_source table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : list of DetectedSource objects
        DetectedSource objects to be inserted into
        the assoc_source database table.

    Returns
    -------
    sources : list of DetectedSource objects
        DetectedSource objects with updated assoc_id
        attribute.
    """
    cur = conn.cursor()
    
    for src in sources:
        cur.execute('''INSERT INTO assoc_source (
            ra, e_ra, dec, e_dec, beam, ndetect) VALUES (
            %s, %s, %s, %s, %s, %s)
            RETURNING id''',
                    (src.ra, src.e_ra, src.dec, src.e_dec,
                     src.beam, src.ndetect))
        src.id = cur.fetchone()[0]
        src.assoc_id = src.id
        cur.execute('''UPDATE detected_source SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()
    
    return sources


def update_matched_assoc(conn, sources):
    """Updates position and size properties of an assoc_source
    entry to the new weigthed average from all detections.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : list of DetectedSource objects
        DetectedSource objects extracted from the
        assoc_source table which have been matched
        to sources detected in the current image.
    """
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE assoc_source SET ra = %s, e_ra = %s, dec = %s,
            e_dec = %s, ndetect = %s WHERE id = %s''',
                    (src.ra, src.e_ra, src.dec, src.e_dec, src.ndetect, src.id))

    conn.commit()
    cur.close()


def update_detected_associd(conn, sources):
    """Updates the assoc_id for sources in the
    detected_source table which have been successfully
    associated with existing VLITE sources in the 
    assoc_source table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : list of DetectedSource objects
        DetectedSource objects detected in the current
        image that have been associated with previous
        VLITE sources in the assoc_source table.
    """
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE detected_source SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()


def update_assoc_nmatches(conn, sources):
    """Updates the number of sky catalog matches to
    a given source in the assoc_source table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : list of DetectedSource objects
        DetectedSource objects with updated nmatches
        attribute after running sky catalog matching.
    """
    cur = conn.cursor()

    for src in sources:
        cur.execute('UPDATE assoc_source SET nmatches = %s WHERE id = %s',
                    (src.nmatches, src.id))

    conn.commit()
    cur.close()


def add_catalog_match(conn, sources):
    """Adds an entry to the catalog_match table for
    every sky catalog source matched to a VLITE source
    in the assoc_source table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : list of CatalogSource objects
        CatalogSource matched to a VLITE DetectedSource.
    """
    cur = conn.cursor()

    for src in sources:
        cur.execute('''INSERT INTO catalog_match (
            catalog_id, src_id, assoc_id, min_deruiter) VALUES (
            %s, %s, %s, %s) ON CONFLICT (
            catalog_id, src_id, assoc_id) DO NOTHING''',
                    (src.catalog_id, src.id, src.assoc_id, src.min_deruiter))

    conn.commit()
    cur.close()


def check_catalog_match(conn, asrc_id, catalog):
    """Checks if a source in the assoc_source table already
    has a match to a source in the specified sky catalog
    in the catalog_match table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    asrc_id : integer
        Row id of the source in the assoc_soure table to
        match to the assoc_id in the vlite_unique table.
    catalog : str
        Name of the sky catalog.

    Returns
    -------
    rowid : integer
        The row id of the entry if it exists. Returns
        ``None`` otherwise.
    """
    cur = conn.cursor()

    cur.execute('''SELECT id FROM catalog_match
        WHERE assoc_id = %s AND catalog_id = (
          SELECT id FROM skycat.catalogs WHERE name = %s)''',
                (asrc_id, catalog))
    rowid = cur.fetchone()

    conn.commit()
    cur.close()

    return rowid


def check_vlite_unique(conn, asrc_id):
    """Checks if a given source from the assoc_source
    table is already in the vlite_unique table. This
    is so that sources don't get added twice to the
    vlite_unique table (once when the nmatches = 0 source
    is pulled from assoc_source table and again if no
    sky catalog match is found) when updating the catalog
    matching results by adding new sky catalogs without
    re-doing the previous catalog matching results.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    asrc_id : integer
        Row id of the source in the assoc_soure table to
        match to the assoc_id in the vlite_unique table.

    Returns
    -------
    existing : list
        Returns the id, image_id, and detected columns of
        the entry with the given assoc_id. Otherwise, 
        returns ``None``.
    """
    cur = conn.cursor()

    cur.execute('''SELECT id, image_id, detected FROM vlite_unique
        WHERE assoc_id = %s''', (asrc_id, ))
    existing = cur.fetchall()

    conn.commit()
    cur.close()
   
    return existing


def add_vlite_unique(conn, src, image_id, update=False):
    """Adds an entry in the vlite_unique table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    src : DetectedSource object
        DetectedSource object with attribute nmatches = 0.
    image_id : integer
        Value of the DetectedSource object image_id attribute.
    update : boolean, optional
        If ``True``, the detected column is updated for
        the existing row with the specified image_id and
        assoc_id. Otherwise, a new row is added.
        Default is ``False``.
    """
    cur = conn.cursor()

    if update:
        cur.execute('''UPDATE vlite_unique SET detected = %s
            WHERE image_id = %s AND assoc_id = %s''',
                    (src.detected, image_id, src.id))
    else:
        cur.execute('''INSERT INTO vlite_unique (
            image_id, assoc_id, detected) VALUES (%s, %s, %s)''',
                    (image_id, src.id, src.detected))

    conn.commit()
    cur.close()

    
def get_image_sources(conn, image_id):
    """Returns a list of sources belonging to a
    particular image from the detected_source table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    image_id : integer
        Value of the DetectedSource object image_id attribute.

    Returns
    -------
    detected_sources : List of DetectedSource objects
        Sources pulled from the detected_source table
        translated from row dictionary objects to
        DetectedSource objects.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    cur.execute('SELECT * FROM detected_source WHERE image_id = %s',
                (image_id, ))
    rows = cur.fetchall()

    conn.commit()
    cur.close()

    detected_sources = []
    for row in rows:
        detected_sources.append(dbclasses.DetectedSource())
        dbclasses.dict2attr(detected_sources[-1], row)

    return detected_sources


def get_associated(conn, sources):
    """Returns a list of sources belonging to a
    particular image from the assoc_source table.
    
    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : List of DetectedSource objects
        DetectedSource objects pulled from the
        detected_source table based on image_id.

    Returns
    -------
    assoc_sources : List of DetectedSource objects
        Sources pulled from the assoc_source table
        based on matching assoc_id and translated
        from row dictionary objects to 
        DetectedSource objects.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    rows = []
    for src in sources:
        cur.execute('SELECT * FROM assoc_source WHERE id = %s',
                    (src.assoc_id, ))
        rows.append(cur.fetchone())

    conn.commit()
    cur.close()

    assoc_sources = []
    for row in rows:
        if not row:
            continue
        assoc_sources.append(dbclasses.DetectedSource())
        dbclasses.dict2attr(assoc_sources[-1], row)

    return assoc_sources


def delete_matches(conn, sources, image_id):
    """Deletes all previous sky catalog cross-matching
    results for a given set of sources. This function
    is called when the config. file option "redo match"
    is ``True``.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : List of DetectedSource objects
        DetectedSource objects belonging to a
        particular image pulled from the
        assoc_source table based on matching assoc_id.
    image_id : integer
        Image id to match to the image_id in the
        vlite_unique table.

    Returns
    -------
    sources : List of DetectedSource objects
        DetectedSource objects with their nmatches
        attribute re-initialized to ``None``.
    """
    print('\nRemoving previous sky catalog matching results '
          'for {} sources.'.format(len(sources)))

    cur = conn.cursor()

    for src in sources:
        src.nmatches = None
        cur.execute('UPDATE assoc_source SET nmatches = %s WHERE id = %s',
                    (src.nmatches, src.id))
        cur.execute('DELETE FROM catalog_match WHERE assoc_id = %s',
                    (src.id, ))
        cur.execute('''DELETE FROM vlite_unique WHERE image_id = %s
            AND assoc_id = %s''', (image_id, src.id))

    conn.commit()
    cur.close()

    return sources


def update_stage(conn, imobj):
    """Updates the stage column in the image table.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Image object used to set table column values.
    """
    cur = conn.cursor()

    cur.execute('UPDATE image SET stage = %s WHERE id = %s',
                (imobj.stage, imobj.id))

    conn.commit()
    cur.close()


def pybdsf_fail(conn, imobj):
    """Updates the image table error_id column to indicate
    a failure to process by PyBDSF: error_id = 3.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    imobj : database.dbclasses.Image instance
        Image object used to set table column values.
    """
    imobj.error_id = 3
    
    cur = conn.cursor()

    # Update Image table error_id code
    cur.execute('''UPDATE image SET error_id = %s WHERE id = %s''',
                (imobj.error_id, imobj.id))

    conn.commit()
    cur.close()
