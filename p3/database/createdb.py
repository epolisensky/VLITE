"""Creates database schema."""


import sqlite3


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
        DROP TABLE IF EXISTS Island;
        DROP TABLE IF EXISTS Source;
        ''')
        print('\nCreating new tables...')
        cur.executescript('''
        CREATE TABLE Image (
            id INTEGER NOT NULL UNIQUE,
            filename TEXT UNIQUE,
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
            mjdtime NUMERIC,
            tau_time NUMERIC,
            duration NUMERIC,
            nsrc INTEGER,
            rms_box TEXT,
            error_id INTEGER,
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
