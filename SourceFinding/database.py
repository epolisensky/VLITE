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

All database flux/brightness units are received as Jy (or Jy/beam) 
and converted to mJy. All database angular size units are received
as degrees and converted to arcseconds.

Post-Processing Pipeline (P3) Stage 2"""


import sqlite3
import re
import sys
import os
from astropy.io import fits
import pybdsf_source


class ImageTable(object):
    """Creates a database Image table object that can be easily 
    passed to SQL table insertion functions. Object attributes
    correspond to the Image table column values."""
    def __init__(self):
        self.name = None
        self.imsize = None
        self.obs_ra = None
        self.obs_dec = None
        self.pixel_scale = None
        self.obj = None
        self.obs_date = None
        self.map_date = None
        self.freq = None
        self.bmaj = None
        self.bmin = None
        self.bpa = None
        self.noise = None
        self.peak = None
        self.config = None
        self.nvis = None
        self.tau_time = None
        self.duration = None
        self.nsrc = None
        self.rms_box = None


    def header_extract(self, hdr):
        """This method does all the heavy lifting of extracting and
        assigning the desired object attributes (table column values)
        from the fits image header and assigns 'None' if the keyword
        is missing in the header."""
        try:
            naxis1 = hdr['NAXIS1']
            naxis2 = hdr['NAXIS2']
            self.imsize = str((naxis1, naxis2)) # pixels
        except KeyError:
            self.imsize = None
        try:
            self.obs_ra = hdr['OBSRA'] # deg
        except KeyError:
            self.obs_ra = None
        try:
            self.obs_dec = hdr['OBSDEC'] # deg
        except KeyError:
            self.obs_dec = None
        try:
            self.pixel_scale = abs(hdr['CDELT1']) * 3600. # arcsec/pixel
        except KeyError:
            try:
                self.pixel_scale = abs(hdr['CDELT2']) * 3600. # arcsec/pixel
            except KeyError:
                self.pixel_scale = None
        try:
            self.obj = hdr['OBJECT']
        except KeyError:
            self.obj = None
        try:
            self.obs_date = hdr['DATE-OBS']
        except KeyError:
            self.obs_date = None
        try:
            self.map_date = hdr['DATE-MAP']
        except KeyError:
            self.map_date = None
        try:
            self.freq = hdr['RESTFREQ'] / 10**6. # MHz
        except KeyError:
            try:
                if hdr['CTYPE3'] == 'FREQ':
                    self.freq = hdr['CRVAL3'] / 10**6. # MHz
                else:
                    self.freq = hdr['CRVAL4'] / 10**6. # MHz
            except KeyError:
                self.freq = None
        try:
            self.bmaj = hdr['BMAJ'] * 3600. # arcsec
            self.bmin = hdr['BMIN'] * 3600. # arcsec
            self.bpa = hdr['BPA'] # deg
        except KeyError:
            try:
                self.bmaj = hdr['CLEANBMJ'] * 3600. # arcsec
                self.bmin = hdr['CLEANBMN'] * 3600. # arcsec
                self.bpa = hdr['CLEANBPA'] # deg
            except KeyError:
                try:
                    # Search for beam params in AIPS history
                    hl = list(hdr['HISTORY'])
                    for line in hl:
                        x = re.findall('BMAJ=\s+([0-9]\S+)', line)
                        y = re.findall('BMIN=\s+([0-9]\S+)', line)
                        z = re.findall('BPA=\s+([0-9]\S+)', line)
                        if len(x) > 0:
                            self.bmaj = float(x[0]) * 3600. # arcsec
                        if len(y) > 0:
                            self.bmin = float(y[0]) * 3600. # arcsec
                        if len(z) > 0:
                            self.bpa = float(z[0]) # deg
                except KeyError:
                    self.bmaj = None
                    self.bmin = None
                    self.bpa = None
        try:
            self.noise = hdr['ACTNOISE'] * 1000. # mJy/beam
        except KeyError:
            self.noise = None
        try:
            self.peak = hdr['PEAK'] * 1000. # mJy/beam
        except KeyError:
            try:
                self.peak = hdr['DATAMAX'] * 1000. # mJy/beam
            except KeyError:
                self.peak = None
        try:
            self.config = hdr['CONFIG']
        except KeyError:
            self.config = None
        try:
            self.nvis = hdr['NVIS']
        except KeyError:
            self.nvis = None
        try:
            self.tau_time = hdr['TAU_TIME'] # sec
        except KeyError:
            self.tau_time = None
        try:
            self.duration = hdr['DURATION'] # sec
        except KeyError:
            self.duration = None


def f(num):
    """Formats float data for database tables."""
    return '{:.6f}'.format(num)


def log_extract(path):
    """Read PyBDSF log output file and extract the last
    written rms_box parameter and residual standard 
    deviation. Takes as input the directory path to the
    log file(s)."""
    with open(path) as f:
        log = f.read()
    x = re.findall('.*rms_box.*:\s(\(.*\))', log)
    rms_box = x[-1]
    y = re.findall('.*std. dev:\s(.*)\s\(', log)
    gresid_std = float(y[-1]) * 1000. # mJy/beam
    return rms_box, gresid_std


def create_db(dbname):
    """Creates a new database by dropping tables if they exist.
    USE WITH CAUTION! DELETING THE TABLES CANNOT BE UNDONE."""
    cont = raw_input(('\nWARNING: Any existing tables in this database '
                         'will be deleted. Continue? '))
    if cont == 'y' or cont == 'yes':
        conn = sqlite3.connect(dbname)
        cur = conn.cursor()
        print('\nDropping tables...')
        cur.executescript('''
        PRAGMA foreign_keys = "1";
        DROP TABLE IF EXISTS Image;
        DROP TABLE IF EXISTS Island;
        DROP TABLE IF EXISTS Source;
        ''')
        print('\nCreating new tables...')
        cur.executescript('''
        CREATE TABLE Image (
            id INTEGER NOT NULL UNIQUE,
            name TEXT UNIQUE,
            imsize TEXT,
            obs_ra REAL,
            obs_dec REAL,
            pixel_scale REAL,
            object TEXT,
            obs_date NUMERIC,
            map_date NUMERIC,
            freq REAL,
            bmaj REAL,
            bmin REAL,
            bpa REAL,
            noise REAL,
            peak REAL,
            config TEXT,
            nvis INTEGER,
            tau_time NUMERIC,
            duration NUMERIC,
            nsrc INTEGER,
            rms_box TEXT,
            PRIMARY KEY(id)
        );

        CREATE TABLE Island (
            isl_id INTEGER NOT NULL,
            image_id INTEGER,
            total_flux REAL,
            e_total_flux REAL,
            rms REAL,
            mean REAL,
            resid_rms REAL,
            resid_mean REAL,
            PRIMARY KEY(isl_id, image_id),
            FOREIGN KEY(image_id) REFERENCES Image(id) ON DELETE CASCADE
        );

        CREATE TABLE Source (
            src_id INTEGER NOT NULL,
            isl_id INTEGER,
            image_id INTEGER,
            ra REAL,
            e_ra REAL,
            dec REAL,
            e_dec REAL,
            total_flux REAL,
            e_total_flux REAL,
            peak_flux REAL,
            e_peak_flux REAL,
            ra_max REAL,
            e_ra_max REAL,
            dec_max REAL,
            e_dec_max REAL,
            maj REAL,
            e_maj REAL,
            min REAL,
            e_min REAL,
            pa REAL,
            e_pa REAL,
            dc_maj REAL,
            e_dc_maj REAL,
            dc_min REAL,
            e_dc_min REAL,
            dc_pa REAL,
            e_dc_pa REAL,
            code TEXT,
            catalog_id INTEGER,
            match_id INTEGER,
            min_deRuiter REAL,
            PRIMARY KEY(src_id, isl_id, image_id),
            FOREIGN KEY(isl_id, image_id) REFERENCES Island(isl_id, image_id) 
              ON DELETE CASCADE
        );
        ''')

        conn.commit()
        cur.close()

    else:
        print("\nNo new database created.")


def filldb_tables(dbname, it, sources):
    """Method to fill existing database tables. Takes as
    input the database name, an ImageTable object, and a 
    list of DetectedSource objects. Attributes of these 
    objects correspond to the table columns. This method 
    first checks if the image name (including the full 
    directory path) already exists. If it does, it UPDATES 
    the Image table entry and DELETES corresponding Island 
    and Source table entries then INSERTS the new values."""
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()
    cur.execute('PRAGMA foreign_keys = "1"')
      
    # Check if image has already been run & insert or update accordingly
    cur.execute('SELECT id FROM Image WHERE name = ?', (it.name, ))
    exists = cur.fetchone()
    if exists is None:
        print("\nAdding {} to database\n".format(it.name))
        cur.execute('''INSERT INTO Image (
            name, imsize, obs_ra, obs_dec, pixel_scale, object, obs_date, 
            map_date, freq, bmaj, bmin, bpa, noise, peak, config, nvis,
            tau_time, duration, nsrc, rms_box) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 
            ?, ?, ?, ?, ?)''', (it.name, it.imsize, it.obs_ra, it.obs_dec,
                                it.pixel_scale, it.obj, it.obs_date,
                                it.map_date, it.freq, it.bmaj, it.bmin, it.bpa,
                                it.noise, it.peak, it.config, it.nvis,
                                it.tau_time, it.duration, it.nsrc, it.rms_box))
        # Get new image id
        cur.execute('SELECT id FROM Image WHERE name = ?', (it.name, ))
        img_id = cur.fetchone()[0]
    else:
        print("\nUpdating existing entries for {}\n".format(it.name))
        img_id = exists[0]
        cur.execute('''UPDATE Image SET name = ?, imsize = ?, obs_ra = ?,
            obs_dec = ?, pixel_scale = ?, object = ?, obs_date = ?, 
            map_date = ?, freq = ?, bmaj = ?, bmin = ?, bpa = ?, noise = ?, 
            peak = ?, config = ?, nvis = ?, tau_time = ?, duration = ?, 
            nsrc = ?, rms_box = ? WHERE id = ?''',
                    (it.name, it.imsize, it.obs_ra, it.obs_dec,
                     it.pixel_scale, it.obj, it.obs_date,
                     it.map_date, it.freq, it.bmaj, it.bmin, it.bpa,
                     it.noise, it.peak, it.config, it.nvis, it.tau_time,
                     it.duration, it.nsrc, it.rms_box, img_id))
        # Delete corresponding Island & Source table entries
        cur.execute('DELETE FROM Island WHERE image_id = ?', (img_id, ))

    for src in sources:
        # Insert values into Island table
        cur.execute('''INSERT OR IGNORE INTO Island (
            isl_id, image_id, total_flux, e_total_flux, 
            rms, mean, resid_rms, resid_mean) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)''',
            (src.isl_id, img_id, f(src.total_flux_isl), f(src.total_flux_islE),
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
            (src.src_id, src.isl_id, img_id, f(src.ra), f(src.e_ra),
             f(src.dec), f(src.e_dec), f(src.total_flux), f(src.e_total_flux),
             f(src.peak_flux), f(src.e_peak_flux), f(src.ra_max),
             f(src.e_ra_max), f(src.dec_max), f(src.e_dec_max), f(src.maj),
             f(src.e_maj), f(src.min), f(src.e_min), f(src.pa), f(src.e_pa),
             f(src.dc_maj), f(src.e_dc_maj), f(src.dc_min), f(src.e_dc_min),
             f(src.dc_pa), f(src.e_dc_pa), src.code, None, None, None))

    conn.commit()
    cur.close()


def pybdsf_to_db(dbname, out):
    """Initializes necessary inputs and then calls 
    filldb_tables() when the object on hand is the output
    from the PyBDSF source finding task process_image()."""
    # Initialize ImageTable object
    it = ImageTable()
    it.name = out.filename
    it.nsrc = out.nsrc
    it.rms_box = str(out.rms_box)
    it.header_extract(out.header)
    if it.imsize is None:
        it.imsize = str(out._original_shape) # pixels
    if it.freq is None:
        it.freq = out.frequency / 10**6. # MHz
    if it.bmaj is None:
        it.bmaj = out.beam[0] # deg
    if it.bmin is None:
        it.bmin = out.beam[1] # deg
    if it.bpa is None:
        it.bpa = out.beam[2] # deg

    # Create new DetectedSource objects from PyBDSF output objects
    newsrcs = []
    for oldsrc in out.sources:
        newsrcs.append(pybdsf_source.DetectedSource())
        newsrcs[-1].cast(oldsrc)

    # Pass the objects to table writing function
    filldb_tables(dbname, it, newsrcs)

    # Return ImageTable & DetectedSource objects for use in cross-matching
    return it, newsrcs
    

def cat_to_db(dbname, cat):
    """Initializes necessary inputs and then calls
    filldb_tables() when the object on hand is a pre-existing
    ascii source catalog. The catalog is read using the 
    pybdsf_source.py task read_catalog, which creates a list of 
    DetectedSource objects. The PyBDSF output log is also accessed
    to input values into the Image table."""
    print('\nReading catalog {}'.format(re.findall('\S+/(\S+)', cat)[0]))
    sources = pybdsf_source.read_catalog(cat)

    # Initialize ImageTable object
    it = ImageTable()
    it.nsrc = len(sources)

    # Find image & log files
    '''Assumes there is a directory called 'catalogs' inside the 
    directory containing the fits images:'''
    path = re.findall('(\S+/)catalogs', cat)[0]
    imname = re.findall('\S+/(\S+).pybdsm', cat)[0] + '.fits'
    it.name = os.path.join(path, imname)
    data, header = fits.getdata(it.name, header=True)
    it.header_extract(header)
    logname = imname + '.pybdsf.log'
    logdir = 'logs/' + logname
    logpath = os.path.join(path, logdir)
    it.rms_box, it.noise = log_extract(logpath)


    # Pass the objects to table writing function
    filldb_tables(dbname, it, sources)

    # Return ImageTable & DetectedSource objects for use in cross-matching
    return it, sources
