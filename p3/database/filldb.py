"""This module contains methods to prepare and insert
data into the ``PostgreSQL`` database tables.

"""
import numpy as np
import re
from astropy.io import fits
import pybdsfcat
import dbclasses


def f(num):
    """Formats float data for database tables."""
    return '{:.6f}'.format(num)


def initImage(impath):
    """Initializes an object of the Image class and sets values 
    for its attributes from the fits file header."""
    img = dbclasses.Image(impath)
    data, header = img.read()
    img.header_attrs(header) # Use header info to set attributes

    return img


def statusCheck(conn, impath):
    """Returns the id and highest completed stage from the 
    database image table or ``None`` if the image filename is 
    not in the database."""
    cur = conn.cursor()
    cur.execute('SELECT id, stage FROM Image WHERE filename = %s', (impath, ))
    status = cur.fetchone()
    cur.close()

    return status


def addImage(conn, impath, status, delete=False):
    """Insert or update rows in database image table.
    
    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL database connection object.
    impath : str
        Directory path to the fits image file.
    status : tuple
        Contains the id and stage for the image's row
        entry in the database image table if it exists.
        Otherwise, status is ``None``.
    delete : bool
        If ``True``, rows in the raw_island database table
        with the appropriate image_id will be deleted,
        cascading to the raw_source table and triggering 
        updates on the assoc_source table after an update
        on the image table. Rows with the same image_id
        are also deleted from the null_detections table.
        Default value is ``False``.

    Returns
    -------
    img : database.dbclasses.Image instance
        Initialized Image object which was used to
        insert values into the database image table.    
    """
    # Initialize Image object
    img = initImage(impath)    
    
    cur = conn.cursor()
    
    # Add new image to DB
    if status is None:
        print('\nAdding {} to database'.format(img.filename))
        sql = '''INSERT INTO image (
            filename, imsize, obs_ra, obs_dec, pixel_scale, object, obs_date, 
            map_date, obs_freq, primary_freq, bmaj, bmin, bpa, noise, peak, 
            config, nvis, mjdtime, tau_time, duration, nsrc, rms_box, 
            error_id, stage) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id'''
        vals = (img.filename, img.imsize, img.obs_ra, img.obs_dec,
                img.pixel_scale, img.obj, img.obs_date, img.map_date,
                img.obs_freq, img.pri_freq, img.bmaj, img.bmin, img.bpa,
                img.noise, img.peak, img.config, img.nvis, img.mjdtime,
                img.tau_time, img.duration, img.nsrc, img.rms_box,
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
            nvis = %s, mjdtime = %s, tau_time = %s, duration = %s, nsrc = %s, 
            rms_box = %s, error_id = %s, stage = %s
            WHERE id = %s'''
        vals = (img.filename, img.imsize, img.obs_ra, img.obs_dec,
                img.pixel_scale, img.obj, img.obs_date, img.map_date,
                img.obs_freq, img.pri_freq, img.bmaj, img.bmin, img.bpa,
                img.noise, img.peak, img.config, img.nvis, img.mjdtime,
                img.tau_time, img.duration, img.nsrc, img.rms_box,
                img.error_id, img.stage, img.id)
        cur.execute(sql, vals)
        if delete:
            # Delete corresponding raw_island & raw_source table entries
            print('\nRemoving previous sources...')
            cur.execute('DELETE FROM raw_island WHERE image_id = %s', (
                img.id, ))
            cur.execute('DELETE FROM null_detections WHERE image_id = %s', (
                img.id, ))

    conn.commit()
    cur.close()

    return img


def addSources(conn, img, sources):
    """Inserts DetectedSource objects into database
    raw_island and raw_source tables. The Image object
    and image table are also updated with some results
    from the source finding.

    """
    cur = conn.cursor()

    # Update Image table -- overwrites all possible updated columns
    sql = '''UPDATE Image SET imsize = %s, obs_freq = %s, bmaj = %s, 
        bmin = %s, bpa = %s, noise = %s, nsrc = %s, rms_box = %s, 
        error_id = %s, stage = %s WHERE id = %s'''
    vals = (img.imsize, img.obs_freq, img.bmaj, img.bmin, img.bpa, img.noise,
            img.nsrc, img.rms_box, img.error_id, img.stage, img.id)
    cur.execute(sql, vals)
    
    # Insert DetectedSources into raw_source and raw_island tables
    print('\nAdding detected sources to database.')
    for src in sources:
        src.image_id = img.id
        # Insert values into raw_island table
        sql = '''INSERT INTO raw_island (
            isl_id, image_id, total_flux, e_total_flux, 
            rms, mean, resid_rms, resid_mean) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (isl_id, image_id)
            DO NOTHING'''
        vals = (src.isl_id, src.image_id, f(src.total_flux_isl),
                f(src.total_flux_islE), f(src.rms_isl), f(src.mean_isl),
                f(src.resid_rms), f(src.resid_mean))
        cur.execute(sql, vals)

        # Insert values into raw_source table
        sql = '''INSERT INTO raw_source (
            src_id, isl_id, image_id, ra, e_ra, dec, e_dec,
            total_flux, e_total_flux, peak_flux, e_peak_flux, 
            ra_max, e_ra_max, dec_max, e_dec_max, maj, e_maj, 
            min, e_min, pa, e_pa, dc_maj, e_dc_maj, dc_min, e_dc_min,
            dc_pa, e_dc_pa, code, assoc_id) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'''
        vals = (src.src_id, src.isl_id, src.image_id, f(src.ra), f(src.e_ra),
                f(src.dec), f(src.e_dec), f(src.total_flux),
                f(src.e_total_flux), f(src.peak_flux), f(src.e_peak_flux),
                f(src.ra_max), f(src.e_ra_max), f(src.dec_max),
                f(src.e_dec_max), f(src.maj), f(src.e_maj), f(src.min),
                f(src.e_min), f(src.pa), f(src.e_pa), f(src.dc_maj),
                f(src.e_dc_maj), f(src.dc_min), f(src.e_dc_min), f(src.dc_pa),
                f(src.e_dc_pa), src.code, src.assoc_id)
        cur.execute(sql, vals)

    conn.commit()
    cur.close()


def addAssoc(dbname, sources):
    print('\nAdding {} newly detected sources to database AssocSource '
          'table'.format(len(sources)))
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()
    
    for src in sources:
        cur.execute('''INSERT INTO AssocSource (
            ra, e_ra, dec, e_dec, maj, e_maj, min, e_min, pa, e_pa, beam,
            ndetect, nopp)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)''',
                    (f(src.ra), f(src.e_ra), f(src.dec),
                     f(src.e_dec), f(src.maj), f(src.e_maj), f(src.min),
                     f(src.e_min), f(src.pa), f(src.e_pa), src.beam,
                     src.ndetect, src.nopp))
        rowid = cur.lastrowid
        cur.execute('''UPDATE RawSource SET assoc_id = %s WHERE src_id = %s AND
            image_id = %s''', (rowid, src.src_id, src.image_id))

    conn.commit()
    cur.close()


def updateMatchedAssoc(dbname, sources):
    print('\nUpdating AssocSource parameters for {} rows'.format(len(sources)))
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE AssocSource SET ra = %s, e_ra = %s, dec = %s,
            e_dec = %s, maj = %s, e_maj = %s, min = %s, e_min = %s, pa =%s,
            e_pa = %s, ndetect = %s WHERE id = %s''',
                    (f(src.ra), f(src.e_ra), f(src.dec), f(src.e_dec),
                     f(src.maj), f(src.e_maj), f(src.min), f(src.e_min),
                     f(src.pa), f(src.e_pa), src.ndetect, src.id))

    conn.commit()
    cur.close()


def updateNullAssoc(dbname, sources, imgid):
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()

    for src in sources:
        cur.execute('UPDATE AssocSource SET nopp = %s WHERE id = %s',
                    (src.nopp, src.id))
        cur.execute('''INSERT INTO NullDetections (assoc_id, image_id)
            VALUES (%s, %s)''', (src.id, imgid))

    conn.commit()
    cur.close()


def updateRawAssocid(dbname, sources):
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE RawSource SET assoc_id = %s WHERE src_id = %s AND
            image_id = %s''', (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()


def pybdsf_fail(dbname, img):
    """Updates the Image table error_id column to indicate
    a failure to process by PyBDSF: error_id = 3."""
    img.error_id = 3
    
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()

    # Update Image table error_id code
    cur.execute('''UPDATE Image SET error_id = %s WHERE filename = %s''',
                (img.error_id, img.filename))

    conn.commit()
    cur.close()

