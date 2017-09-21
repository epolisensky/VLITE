"""Creates database schema."""


import sqlite3


def makeError(cursor):
    reason_dict = {'short duration' : 1,
                   'close to bright radio source' : 2,
                   'PyBDSF failed to process' : 3,
                   'poor quality' : 4}

    for key in reason_dict.keys():
        cursor.execute('INSERT INTO Error (id, reason) VALUES (?, ?)',
                       (reason_dict[key], key))


def create(dbname, safe):
    """Creates a new database by dropping tables if they exist.
    USE WITH CAUTION! DELETING THE TABLES CANNOT BE UNDONE."""
    if not safe:
        cont = raw_input(('\nWARNING: Any existing tables in this database '
                          'will be deleted. Continue? '))
    else:
        cont = 'yes'
    if cont == 'y' or cont == 'yes':
        print('\nMaking new database...')
        conn = sqlite3.connect(dbname)
        cur = conn.cursor()
        print('\nDropping tables...')
        cur.executescript('''
        PRAGMA foreign_keys = "1";
        DROP TABLE IF EXISTS Image;
        DROP TABLE IF EXISTS rawIsland;
        DROP TABLE IF EXISTS rawSource;
        DROP TABLE IF EXISTS Error;
        DROP TABLE IF EXISTS AssocSource
        ''')
        print('\nCreating new tables...')
        cur.executescript('''
        CREATE TABLE Image (
            id INTEGER NOT NULL UNIQUE,
            filename TEXT UNIQUE,
            imsize BLOB,
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
            mjdtime INTEGER,
            tau_time INTEGER,
            duration INTEGER,
            nsrc INTEGER,
            rms_box BLOB,
            error_id INTEGER,
            PRIMARY KEY(id),
            FOREIGN KEY(error_id) REFERENCES Error(id)
        );

        CREATE TABLE Error (
            id INTEGER NOT NULL,
            reason TEXT,
            PRIMARY KEY(id)
        );

        CREATE TABLE rawIsland (
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

        CREATE TABLE rawSource (
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
            assoc_id INTEGER,
            PRIMARY KEY(src_id, image_id),
            FOREIGN KEY(isl_id, image_id) REFERENCES 
              rawIsland(isl_id, image_id) ON DELETE CASCADE,
            FOREIGN KEY(assoc_id) REFERENCES AssocSource(id)
        );

        CREATE TABLE AssocSource (
            id INTEGER NOT NULL,
            ra REAL,
            e_ra REAL,
            dec REAL,
            e_dec REAL,
            maj REAL,
            e_maj REAL,
            min REAL,
            e_min REAL,
            pa REAL,
            e_pa REAL,
            num_detect INTEGER,
            num_null INTEGER,
            catalog_id INTEGER,
            match_id INTEGER,
            min_deRuiter REAL,
            PRIMARY KEY(id)
        );
        ''')

        makeError(cur)

        conn.commit()
        cur.close()

    else:
        print("\nNo new database created.")
