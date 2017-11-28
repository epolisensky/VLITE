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


def init_image(impath):
    """Initializes an object of the Image class and sets values 
    for its attributes from the fits file header."""
    img = dbclasses.Image(impath)
    data, header = img.read()
    img.header_attrs(header) # Use header info to set attributes

    return img


def status_check(conn, impath):
    """Returns the id and highest completed stage from the 
    database image table or ``None`` if the image filename is 
    not in the database."""
    cur = conn.cursor()
    cur.execute('SELECT id, stage FROM Image WHERE filename = %s', (impath, ))
    status = cur.fetchone()
    cur.close()

    return status


def add_image(conn, impath, status, delete=False):
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
        If ``True``, rows in the detected_island database table
        with the appropriate image_id will be deleted,
        cascading to the detected_source table and triggering 
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
    img = init_image(impath)    
    
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

    """
    cur = conn.cursor()

    # Update image table -- overwrites all possible updated columns
    sql = '''UPDATE image SET imsize = %s, obs_freq = %s, bmaj = %s, 
        bmin = %s, bpa = %s, noise = %s, nsrc = %s, rms_box = %s, 
        error_id = %s, stage = %s WHERE id = %s'''
    vals = (img.imsize, img.obs_freq, img.bmaj, img.bmin, img.bpa, img.noise,
            img.nsrc, img.rms_box, img.error_id, img.stage, img.id)
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
        vals = (src.isl_id, src.image_id, f(src.total_flux_isl),
                f(src.total_flux_islE), f(src.rms_isl), f(src.mean_isl),
                f(src.resid_rms), f(src.resid_mean))
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


def add_assoc(conn, sources):
    cur = conn.cursor()
    
    for src in sources:
        cur.execute('''INSERT INTO assoc_source (
            ra, e_ra, dec, e_dec, maj, e_maj, min, e_min, pa, e_pa, beam,
            ndetect) VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id''',
                    (f(src.ra), f(src.e_ra), f(src.dec), f(src.e_dec),
                     f(src.maj), f(src.e_maj), f(src.min), f(src.e_min),
                     f(src.pa), f(src.e_pa), src.beam, src.ndetect))
        src.id = cur.fetchone()[0]
        src.assoc_id = src.id
        cur.execute('''UPDATE detected_source SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()
    
    return sources


def update_matched_assoc(conn, sources):
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE assoc_source SET ra = %s, e_ra = %s, dec = %s,
            e_dec = %s, maj = %s, e_maj = %s, min = %s, e_min = %s, pa =%s,
            e_pa = %s, beam = %s, ndetect = %s WHERE id = %s''',
                    (f(src.ra), f(src.e_ra), f(src.dec), f(src.e_dec),
                     f(src.maj), f(src.e_maj), f(src.min), f(src.e_min),
                     f(src.pa), f(src.e_pa), src.beam, src.ndetect, src.id))

    conn.commit()
    cur.close()


def update_detected_associd(conn, sources):
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE detected_source SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()


def update_assoc_nmatches(conn, sources):
    cur = conn.cursor()

    for src in sources:
        cur.execute('UPDATE assoc_source SET nmatches = %s WHERE id = %s',
                    (src.nmatches, src.id))

    conn.commit()
    cur.close()


def add_catalog_match(conn, sources):
    cur = conn.cursor()

    for src in sources:
        cur.execute('''INSERT INTO catalog_match (
            catalog_id, src_id, assoc_id, min_deruiter) VALUES (
            %s, %s, %s, %s)''',
                    (src.catalog_id, src.id, src.assoc_id, src.min_deruiter))

    conn.commit()
    cur.close()


def add_vlite_unique(conn, src, image_id):
    cur = conn.cursor()

    cur.execute('''INSERT INTO vlite_unique (
        image_id, assoc_id, detected) VALUES (%s, %s, %s)''',
                (image_id, src.id, src.detected))

    conn.commit()
    cur.close()


def update_stage(conn, imobj):
    cur = conn.cursor()

    cur.execute('UPDATE image SET stage = %s WHERE id = %s',
                (imobj.stage, imobj.id))

    conn.commit()
    cur.close()


def pybdsf_fail(conn, imobj):
    """Updates the image table error_id column to indicate
    a failure to process by PyBDSF: error_id = 3.

    """
    imobj.error_id = 3
    
    cur = conn.cursor()

    # Update Image table error_id code
    cur.execute('''UPDATE image SET error_id = %s WHERE id = %s''',
                (imobj.error_id, imobj.id))

    conn.commit()
    cur.close()

