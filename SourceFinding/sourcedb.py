"""sourcedb.py creates/updates an SQLite database to store
results from source finding on an image using PyBDSF (see
runPyBDSF.py). The database contains three tables: source 
parameters (location, flux, size, etc.) are inserted into
the Source table, properties of the islands in which the
sources reside, such as flux, rms, mean, and residual values
are inserted into the Island table, and metadata characterizing
the images are inserted into the Image table. Source & island
values are extracted from the bdsf.process_image() output object 
or from a pre-existing PyBDSF ascii source catalog. Image values
primarily come from the image header and PyBDSF log output. An
Image table object is created to more easily pass all these values 
to the filldb_tables() method.

All database flux/brightness units are received as Jy (or Jy/beam) 
and converted to mJy. All database angular size units are received
as degrees and converted to arcseconds.

Post-Processing Pipeline (P3) Stage 2"""


import sqlite3
import re
import sys
import os
from astropy.io import fits
import readPyBDSFcat


class ImageTable(object):
    """Creates a database Image table object that can be easily 
    passed to SQL table insertion functions. Object attributes
    correspond to the Image table column values."""
    def __init__(self):
        self.name = ''
        self.imsize = '(-999, -999)' # pixels
        self.obs_ra = -99.9 # deg
        self.obs_dec = -99.9 # deg
        self.pixscale = -99.9 # deg --> gets converted to arcsec
        self.obj = ''
        self.obs_date = 'yyyy-mm-dd' # date format will vary
        self.map_date = 'yyyy-mm-dd' # date format will vary
        self.freq = -99.9 # MHz
        self.bmaj = -99.9 # deg --> gets converted to arcsec
        self.bmin = -99.9 # deg --> gets converted to arcsec
        self.bpa = -99.9 # deg
        self.noise = -99.9 # Jy/beam --> gets converted to mJy/beam
        self.peak = -99.9 # Jy/beam --> gets converted to mJy/beam
        self.config = ''
        self.nvis = -9999
        self.tau_time = -9999 # seconds? integration time?
        self.duration = -9999 # seconds? total TOS?
        self.nsrc = -9999 # number of sources found by PyBDSF
        self.rms_box = '(-999, -999)' # box size used by PyBDSF


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
            self.pixscale = abs(hdr['CDELT1']) # deg
        except KeyError:
            try:
                self.pixscale = abs(hdr['CDELT2']) # deg
            except KeyError:
                self.pixscale = None
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
            self.bmaj = hdr['BMAJ'] # deg
            self.bmin = hdr['BMIN'] # deg
            self.bpa = hdr['BPA'] # deg
        except KeyError:
            try:
                self.bmaj = hdr['CLEANBMJ'] # deg
                self.bmin = hdr['CLEANBMN'] # deg
                self.bpa = hdr['CLEANBPA'] # deg
            except KeyError:
                try:
                    # search for beam params in AIPS history
                    hl = list(hdr['HISTORY'])
                    for line in hl:
                        x = re.findall('BMAJ=\s+([0-9]\S+)', line)
                        y = re.findall('BMIN=\s+([0-9]\S+)', line)
                        z = re.findall('BPA=\s+([0-9]\S+)', line)
                        if len(x) > 0:
                            self.bmaj = float(x[0]) # deg
                        if len(y) > 0:
                            self.bmin = float(y[0]) # deg
                        if len(z) > 0:
                            self.bpa = float(z[0]) # deg
                except KeyError:
                    self.bmaj = None
                    self.bmin = None
                    self.bpa = None
        try:
            self.noise = hdr['ACTNOISE'] # Jy/beam
        except KeyError:
            self.noise = None
        try:
            self.peak = hdr['PEAK'] # Jy/beam
        except KeyError:
            try:
                self.peak = hdr['DATAMAX'] # Jy/beam
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
    gresid_std = float(y[-1]) # Jy/beam
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
    input the database name, an Image table object (from
    classes.py), and a list of source objects. Attributes
    of these objects should correspond to the table columns.
    This method first checks if the image name (including
    the full directory path) already exists. If it does,
    it UPDATES the Image table entry and DELETES corresponding
    Island & Source table entries then INSERTS the new values."""
    conn = sqlite3.connect(dbname)
    cur = conn.cursor()
    cur.execute('PRAGMA foreign_keys = "1"')
      
    # check if image has already been run & insert or update accordingly
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
                                (it.pixscale * 3600.), it.obj, it.obs_date,
                                it.map_date, it.freq, (it.bmaj * 3600.),
                                (it.bmin * 3600.), it.bpa, (it.noise * 1000.),
                                (it.peak * 1000.), it.config, it.nvis,
                                it.tau_time, it.duration, it.nsrc, it.rms_box))
        # get new image id
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
                     (it.pixscale * 3600.), it.obj, it.obs_date, it.map_date,
                     it.freq, (it.bmaj * 3600.), (it.bmin * 3600.), it.bpa,
                     (it.noise * 1000.), (it.peak * 1000.), it.config,
                     it.nvis, it.tau_time, it.duration,
                     it.nsrc, it.rms_box, img_id))
        # delete corresponding Island & Source table entries
        cur.execute('DELETE FROM Island WHERE image_id = ?', (img_id, ))

    for src in sources:
        # Insert values into Island table
        cur.execute('''INSERT OR IGNORE INTO Island (
            isl_id, image_id, total_flux, e_total_flux, 
            rms, mean, resid_rms, resid_mean) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)''',
            (src.island_id, img_id, f(src.total_flux_isl * 1000.),
             f(src.total_flux_islE * 1000.),
             f(src.rms_isl * 1000.), f(src.mean_isl * 1000.),
             f(src.gresid_rms * 1000.), f(src.gresid_mean * 1000.)))

        # Insert values into Source table
        cur.execute('''INSERT INTO Source (
            src_id, isl_id, image_id, ra, e_ra, dec, e_dec,
            total_flux, e_total_flux, peak_flux, e_peak_flux, 
            ra_max, e_ra_max, dec_max, e_dec_max, maj, e_maj, 
            min, e_min, pa, e_pa, dc_maj, e_dc_maj, dc_min, e_dc_min,
            dc_pa, e_dc_pa, code) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 
              ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
            (src.source_id, src.island_id, img_id,
             f(src.posn_sky_centroid[0]), f(src.posn_sky_centroidE[0]),
             f(src.posn_sky_centroid[1]), f(src.posn_sky_centroidE[1]),
             f(src.total_flux * 1000.), f(src.total_fluxE * 1000.),
             f(src.peak_flux_centroid * 1000.),
             f(src.peak_flux_centroidE * 1000.),
             f(src.posn_sky_max[0]), f(src.posn_sky_maxE[0]),
             f(src.posn_sky_max[1]), f(src.posn_sky_maxE[1]),
             f(src.size_sky[0] * 3600.), f(src.size_skyE[0] * 3600.),
             f(src.size_sky[1] * 3600.), f(src.size_skyE[1] * 3600.),
             f(src.size_sky[2]), f(src.size_skyE[2]),
             f(src.deconv_size_sky[0] * 3600.),
             f(src.deconv_size_skyE[0] * 3600.),
             f(src.deconv_size_sky[1] * 3600.),
             f(src.deconv_size_skyE[1] * 3600.),
             f(src.deconv_size_sky[2]), f(src.deconv_size_skyE[2]),
             src.code))

    conn.commit()
    cur.close()


def pybdsf_to_db(dbname, out):
    """Initializes necessary inputs and then calls 
    filldb_tables() when the object on hand is the output
    from the PyBDSF source finding task process_image()."""
    # initialize Image table object
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

    # pass Image table object & list of source objects
    filldb_tables(dbname, it, out.sources)
    

def cat_to_db(dbname, cat):
    """Initializes necessary inputs and then calls
    filldb_tables() when the object on hand is a pre-existing
    ascii source catalog. The catalog is read using the iofuncs.py
    task readPybdsfCatFile, which creates a list of pybdsf_source
    objects. Attributes are edited/added to this object to match
    the attributes of the bdsf.process_image() output object, which
    has the correct formatting for input into the Source & Island 
    tables. In this case, the PyBDSF output log is also accessed
    to input values into the Image table."""
    print('\nReading catalog {}'.format(re.findall('\S+/(\S+)', cat)[0]))
    src_obj_list = readPyBDSFcat.read_catalog(cat)

    # initialize Image table object
    it = ImageTable()
    it.nsrc = len(src_obj_list)

    # Find image file
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

    key_map = {'source_id' : 'srcID',
               'island_id' : 'islID',
               'total_flux' : 'Total',
               'total_fluxE' : 'dTotal',
               'peak_flux_centroid' : 'Peak',
               'peak_flux_centroidE' : 'dPeak',
               'total_flux_isl' : 'islTotal',
               'total_flux_islE' : 'dislTotal',
               'rms_isl' : 'islRMS',
               'mean_isl' : 'islMEAN',
               'gresid_rms' : 'RRMS',
               'gresid_mean' : 'RMEAN',
               'code': 'SCode'}

    for src in src_obj_list:
        # rename attrs
        for key in key_map.keys():
            src.__dict__[key] = src.__dict__.pop(key_map[key])

        # combine old attrs into new ones
        src.__dict__['posn_sky_centroid'] = [src.__dict__['RA'],
                                             src.__dict__['Dec']]
        src.__dict__['posn_sky_centroidE'] = [src.__dict__['dRA'],
                                              src.__dict__['dDec']]
        src.__dict__['posn_sky_max'] = [src.__dict__['RAmax'],
                                        src.__dict__['Decmax']]
        src.__dict__['posn_sky_maxE'] = [src.__dict__['dRAmax'],
                                         src.__dict__['dDecmax']]
        src.__dict__['size_sky'] = [src.__dict__['BMAJ'],
                                    src.__dict__['BMIN'],
                                    src.__dict__['PA']]
        src.__dict__['size_skyE'] = [src.__dict__['dBMAJ'],
                                     src.__dict__['dBMIN'],
                                     src.__dict__['dPA']]
        src.__dict__['deconv_size_sky'] = [src.__dict__['BMAJDC'],
                                           src.__dict__['BMINDC'],
                                           src.__dict__['PADC']]
        src.__dict__['deconv_size_skyE'] = [src.__dict__['dBMAJDC'],
                                            src.__dict__['dBMINDC'],
                                            src.__dict__['dPADC']]


    # pass Image table object & list of source objects
    filldb_tables(dbname, it, src_obj_list)
