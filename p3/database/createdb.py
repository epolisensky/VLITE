"""Functions to create database tables, functions, and triggers."""


def make_error(cursor):
    """Inserts values into the database error table."""
    reason_dict = {'short duration' : 1,
                   'close to bright radio source' : 2,
                   'PyBDSF failed to process' : 3,
                   'poor quality' : 4}

    sql = 'INSERT INTO error (id, reason) VALUES (%s, %s);'
    for key in reason_dict.keys():
        cursor.execute(sql, (reason_dict[key], key))


def create(conn, safe=False):
    """
    Creates new tables and triggers for the connected `PostgreSQL` 
    database by dropping tables if they exist. The current user
    must own the tables or be a superuser in order to drop them.
    USE WITH CAUTION! DELETING THE DATA CANNOT BE UNDONE.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL database connection object.
    safe : bool, optional
        If ``False``, the user will be warned that existing data
        is about to be deleted and prompted to continue. Default
        value is ``False``.
    """
    if not safe:
        cont = raw_input(('\nWARNING: Any existing tables and data in this '
                          'database will be deleted. Are you sure you want '
                          'to continue? '))
    else:
        cont = 'yes'
    if cont == 'y' or cont == 'yes':
        print('\nDropping tables if they exist...')
        cur = conn.cursor()
        sql = (
            '''
            DROP TABLE IF EXISTS catalog_match;
            DROP TABLE IF EXISTS corrected_flux;
            DROP TABLE IF EXISTS detected_source;
            DROP TABLE IF EXISTS detected_island;
            DROP TABLE IF EXISTS image;
            DROP TABLE IF EXISTS assoc_source;
            DROP TABLE IF EXISTS error;
            DROP FUNCTION IF EXISTS update_assoc_func
            ''')
        cur.execute(sql)

        print('\nCreating new tables...')
        sql = (
            '''
            CREATE EXTENSION IF NOT EXISTS q3c;

            CREATE TABLE error (
                id INTEGER NOT NULL,
                reason TEXT,
                PRIMARY KEY (id)
            );

            CREATE TABLE assoc_source (
                id SERIAL NOT NULL,
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
                beam REAL,
                ndetect INTEGER,
                nmatches INTEGER,
                PRIMARY KEY (id)
            )
            WITH (fillfactor=90);

            CREATE TABLE image (
                id SERIAL NOT NULL UNIQUE,
                filename TEXT UNIQUE,
                imsize VARCHAR(14),
                obs_ra REAL,
                obs_dec REAL,
                pixel_scale REAL,
                object TEXT,
                obs_date DATE,
                map_date DATE,
                obs_freq REAL,
                primary_freq REAL,
                bmaj REAL,
                bmin REAL,
                bpa REAL,
                noise REAL,
                peak REAL,
                config TEXT,
                nvis INTEGER,
                mjdtime REAL,
                tau_time REAL,
                duration REAL,
                nsrc INTEGER,
                rms_box VARCHAR(14),
                stage INTEGER,
                error_id INTEGER,
                PRIMARY KEY (id),
                FOREIGN KEY (error_id) 
                  REFERENCES error (id) 
                  ON UPDATE CASCADE
            );

            CREATE TABLE detected_island (
                isl_id INTEGER NOT NULL,
                image_id INTEGER NOT NULL,
                total_flux REAL,
                e_total_flux REAL,
                rms REAL,
                mean REAL,
                resid_rms REAL,
                resid_mean REAL,
                PRIMARY KEY (isl_id, image_id),
                FOREIGN KEY (image_id) 
                    REFERENCES image (id) 
                    ON DELETE CASCADE
            );

            CREATE TABLE detected_source (
                src_id INTEGER NOT NULL,
                isl_id INTEGER NOT NULL,
                image_id INTEGER NOT NULL,
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
                PRIMARY KEY (src_id, image_id),
                FOREIGN KEY (isl_id, image_id)
                  REFERENCES detected_island (isl_id, image_id)
                  ON DELETE CASCADE,
                FOREIGN KEY (assoc_id)
                  REFERENCES assoc_source (id)
            );

            CREATE TABLE corrected_flux (
                src_id INTEGER NOT NULL,
                isl_id INTEGER NOT NULL,
                image_id INTEGER NOT NULL,
                total_flux REAL,
                e_total_flux REAL,
                peak_flux REAL,
                e_peak_flux REAL,
                isl_total_flux REAL,
                isl_e_total_flux REAL,
                isl_rms REAL,
                isl_mean REAL,
                isl_resid_rms REAL,
                isl_resid_mean REAL,
                PRIMARY KEY (src_id, image_id),
                FOREIGN KEY (src_id, image_id)
                  REFERENCES detected_source (src_id, image_id)
                  ON DELETE CASCADE,
                FOREIGN KEY (isl_id, image_id)
                  REFERENCES detected_island (isl_id, image_id)
                  ON DELETE CASCADE
            );

            CREATE TABLE catalog_match (
                id SERIAL NOT NULL,
                catalog_id INTEGER,
                src_id INTEGER,
                assoc_id INTEGER,
                min_deRuiter REAL,
                PRIMARY KEY (id),
                FOREIGN KEY (assoc_id)
                  REFERENCES assoc_source (id)
                  ON DELETE CASCADE
            );

            CREATE INDEX ON detected_source (q3c_ang2ipix(ra, dec))
                WITH (fillfactor = 90);
            CREATE INDEX ON assoc_source (q3c_ang2ipix(ra, dec))
                WITH (fillfactor = 90);
            ''')
        cur.execute(sql)
        
        # Make error table
        make_error(cur)

        conn.commit()

        # Triggers       
        # Re-compute averages or remove assoc_source after
        # detected_source is deleted
        sql = (
            '''
            CREATE OR REPLACE FUNCTION update_assoc_func()
              RETURNS trigger AS
            $$
            BEGIN
              DELETE FROM assoc_source
              WHERE id = OLD.assoc_id AND ndetect = 1;
              UPDATE assoc_source SET ra = (ra*ndetect-OLD.ra)/(ndetect-1),
                e_ra = (e_ra*ndetect-OLD.e_ra)/(ndetect-1),
                dec = (dec*ndetect-OLD.dec)/(ndetect-1),
                e_dec = (e_dec*ndetect-OLD.e_dec)/(ndetect-1),
                maj = (maj*ndetect-OLD.maj)/(ndetect-1),
                e_maj = (e_maj*ndetect-OLD.e_maj)/(ndetect-1),
                min = (min*ndetect-OLD.min)/(ndetect-1),
                e_min = (e_min*ndetect-OLD.e_min)/(ndetect-1),
                pa = (pa*ndetect-OLD.pa)/(ndetect-1),
                e_pa = (e_pa*ndetect-OLD.e_pa)/(ndetect-1),
                ndetect = ndetect - 1
              WHERE id = OLD.assoc_id;
            RETURN NEW;
            END; $$
            LANGUAGE plpgsql;

            CREATE TRIGGER update_assoc
              AFTER DELETE ON detected_source
              FOR EACH ROW
              EXECUTE PROCEDURE update_assoc_func();
            ''')

        cur.execute(sql)
        conn.commit()
        cur.close()

    else:
        print('\nAborting... database {} left unchanged.'.format(dbname))

