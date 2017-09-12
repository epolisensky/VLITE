"""database.py creates/updates an SQLite database to store
results from source finding on an image using PyBDSF (see
runPyBDSF.py). The database contains three tables: source 
parameters (location, flux, size, etc.) are inserted into
the Source table, properties of the islands in which the
sources reside, such as flux, rms, mean, and residual values
are inserted into the Island table, and metadata characterizing
the images are inserted into the Image table. Source & island
values are extracted from the bdsf.process_image() output object 
or from a pre-existing PyBDSF ascii source catalog and stored in
an DetectedSource object. Image values primarily come from the 
image header and PyBDSF log output. An ImageTable object is created 
to more easily pass all these values to the filldb_tables() method.
It's important that the DetectedSource and ImageTable object
attribute names match the database table column names to allow
consistent access to source properties later on when input to
catalog cross-matching functions can accept a variety of inputs 
(i.e. internal objects with attributes or row dictionaries read
in from external database tables).

All database flux/brightness units mJy or mJy/beam. All angular
size units are arcseconds.

Post-Processing Pipeline (P3) Stage 2"""


import sqlite3
import numpy as np
from astropy.io import fits
import pybdsfcat
import dbclasses
#from quality import test


def f(num):
    """Formats float data for database tables."""
    return '{:.6f}'.format(num)


def initImage(impath):
    img = dbclasses.Image(impath)
    data, header = img.read()
    img.header_attrs(header) # Use header info to set attributes

    # Quality check here

    return img


def addImage(dbname, impath):
    # Initialize Image object
    img = initImage(impath)
    
    """Insert or update database Image table."""
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()
    cur.execute('PRAGMA foreign_keys = "1"') 
    
    # Check if image has already been run and insert or update accordingly
    cur.execute('SELECT id FROM Image WHERE filename = ?', (img.filename, ))
    exists = cur.fetchone()
    if exists is None:
        print("\nAdding {} to database.\n".format(img.filename))
        cur.execute('''INSERT INTO Image (
            filename, imsize, obs_ra, obs_dec, pixel_scale, object, obs_date, 
            map_date, freq, bmaj, bmin, bpa, noise, peak, config, nvis,
            mjdtime, tau_time, duration, nsrc, rms_box, error_id) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 
            ?, ?, ?, ?, ?, ?, ?)''',
                    (img.filename, img.imsize, img.obs_ra, img.obs_dec,
                     img.pixel_scale, img.obj, img.obs_date, img.map_date,
                     img.freq, img.bmaj, img.bmin, img.bpa, img.noise,
                     img.peak, img.config, img.nvis, img.mjdtime, img.tau_time,
                     img.duration, img.nsrc, img.rms_box, img.error_id))
    else:
        print("\nUpdating existing entries for {}\n".format(img.filename))
        imgid = exists[0]
        cur.execute('''UPDATE Image SET filename = ?, imsize = ?, obs_ra = ?,
            obs_dec = ?, pixel_scale = ?, object = ?, obs_date = ?, 
            map_date = ?, freq = ?, bmaj = ?, bmin = ?, bpa = ?, noise = ?, 
            peak = ?, config = ?, nvis = ?, mjdtime = ?, tau_time = ?, 
            duration = ?, nsrc = ?, rms_box = ?, error_id = ? WHERE id = ?''',
                    (img.filename, img.imsize, img.obs_ra, img.obs_dec,
                     img.pixel_scale, img.obj, img.obs_date, img.map_date,
                     img.freq, img.bmaj, img.bmin, img.bpa, img.noise,
                     img.peak, img.config, img.nvis, img.mjdtime, img.tau_time,
                     img.duration, img.nsrc, img.rms_box, img.error_id, imgid))
        # Delete corresponding Island & Source table entries
        cur.execute('DELETE FROM Island WHERE image_id = ?', (imgid, ))

    conn.commit()
    cur.close()

    return img


def addSources(dbname, img, sources):
    """Inserts DetectedSource objects into database Island
    and Source tables. The Image object and table is also
    updated with some results from the source finding."""
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()
    cur.execute('PRAGMA foreign_keys = "1"')

    # Update Image table -- overwrites all possible updated columns
    cur.execute('SELECT id FROM Image WHERE filename = ?', (img.filename, ))
    imgid = cur.fetchone()[0]
    cur.execute('''UPDATE Image SET imsize = ?, freq = ?, bmaj = ?, bmin = ?, 
        bpa = ?, noise = ?, nsrc = ?, rms_box = ?, error_id = ? WHERE id = ?''',
                (img.imsize, img.freq, img.bmaj, img.bmin, img.bpa, img.noise,
                 img.nsrc, img.rms_box, img.error_id, imgid))     

    # Insert DetectedSources into Source and Island tables
    print('\nAdding sources to database.')
    for src in sources:
        # Insert values into Island table
        cur.execute('''INSERT OR IGNORE INTO Island (
            isl_id, image_id, total_flux, e_total_flux, 
            rms, mean, resid_rms, resid_mean) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)''',
            (src.isl_id, imgid, f(src.total_flux_isl), f(src.total_flux_islE),
             f(src.rms_isl), f(src.mean_isl), f(src.resid_rms),
             f(src.resid_mean)))

        # Insert values into Source table
        cur.execute('''INSERT INTO Source (
            src_id, isl_id, image_id, ra, e_ra, dec, e_dec,
            total_flux, e_total_flux, peak_flux, e_peak_flux, 
            ra_max, e_ra_max, dec_max, e_dec_max, maj, e_maj, 
            min, e_min, pa, e_pa, dc_maj, e_dc_maj, dc_min, e_dc_min,
            dc_pa, e_dc_pa, code, catalog_id, match_id, min_deRuiter) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 
              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
            (src.src_id, src.isl_id, imgid, f(src.ra), f(src.e_ra),
             f(src.dec), f(src.e_dec), f(src.total_flux), f(src.e_total_flux),
             f(src.peak_flux), f(src.e_peak_flux), f(src.ra_max),
             f(src.e_ra_max), f(src.dec_max), f(src.e_dec_max), f(src.maj),
             f(src.e_maj), f(src.min), f(src.e_min), f(src.pa), f(src.e_pa),
             f(src.dc_maj), f(src.e_dc_maj), f(src.dc_min), f(src.e_dc_min),
             f(src.dc_pa), f(src.e_dc_pa), src.code, None, None, None))

    conn.commit()
    cur.close()


def pipe_to_table(dbname, img, out):
    """Method to translate PyBDSF output within the
    pipeline to DetectedSource objects and update Image 
    object."""
    # Add PyBDSF defined attributes
    img.nsrc = out.nsrc
    img.rms_box = str(out.rms_box)
    # Try updating any missing attributes from header info
    # using PyBDSF's output object
    if img.imsize is None:
        img.imsize = str(out._original_shape) # pixels
    if img.freq is None:
        img.freq = out.frequency / 10**6. # MHz
    if img.bmaj is None:
        img.bmaj = out.beam[0] * 3600. # arcsec
    if img.bmin is None:
        img.bmin = out.beam[1] * 3600. # arcsec
    if img.bpa is None:
        img.bpa = out.beam[2] # deg
    if img.noise is None:
        img.noise = np.std(out.resid_gaus_arr) * 1000. # mJy/beam

    # Translate PyBDSF output source objects to our own
    # DetectedSource objects
    newsrcs = []
    for oldsrc in out.sources:
        newsrcs.append(dbclasses.DetectedSource())
        newsrcs[-1].cast(oldsrc)

    # Second quality check here
        
    # Pass the objects to table writing function
    addSources(dbname, img, newsrcs)

    # Return Image & DetectedSource objects for use in cross-matching
    return img, newsrcs
    

def cat_to_table(dbname, img, pybdsfdir):
    # Set the file path
    prefix = re.findall('.*\/(.*)\.', img.filename)[0]
    srl = prefix + '.pybdsm.srl'
    catalog = os.path.join(pybdsfdir, srl)
    # Read the catalog
    print('\nExtracting sources from {}'.format(srl))
    sources = pybdsfcat.read_catalog(catalog)

    # Count the number of sources
    img.nsrc = len(sources)

    # Extract rms_box size and noise from log file
    img.rms_box, img.noise, raw_rms = img.log_attrs(pybdsfdir)

    # Pass the objects to table writing function
    addSources(dbname, img, sources)

    # Return ImageTable & DetectedSource objects for use in cross-matching
    return img, sources


def pybdsf_fail(dbname, img):
    """Updates the Image table error_id column to indicate
    a failure to process by PyBDSF: error_id = 3."""
    img.error_id = 3
    
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()

    # Update Image table error_id code
    cur.execute('''UPDATE Image SET error_id = ? WHERE filename = ?''',
                (img.error_id, img.filename))

    conn.commit()
    cur.close()

