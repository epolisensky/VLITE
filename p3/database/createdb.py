"""Creates database tables, functions, and triggers."""
import sys


def make_error(cur, params):
    """Inserts values into the database error table."""
    reason_dict = {'image missing necessary header keyword(s)' : 1,
                   'integration time on source < {} s'.
                   format(params['min time on source (s)']) : 2,
                   'image noise < 0 or > {} mJy/beam'.
                   format(params['max noise (mJy/beam)']): 3,
                   'beam axis ratio > {}'.
                   format(params['max beam axis ratio']) : 4,
                   'image center within {} deg of problem source'.
                   format(params['min problem source separation (deg)']): 5,
                   'PyBDSF failed to process' : 6,
                   'zero sources extracted' : 7,
                   'source metric > {}'.
                   format(params['max source metric']) : 8}

    sql = 'INSERT INTO error (id, reason) VALUES (%s, %s);'
    for key, value in sorted(reason_dict.items(), key=lambda x: x[1]):
        cur.execute(sql, (value, key))



def create(conn, params, safe=False):
    """Creates new tables and triggers for the connected `PostgreSQL` 
    database by dropping tables if they exist. The current user
    must own the tables or be a superuser in order to drop them.
    USE WITH CAUTION! DELETING THE DATA CANNOT BE UNDONE.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    params : dict
        Dictionary of quality check requirements from the run
        configuration file.
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
            DROP TABLE IF EXISTS new_vu;
            DROP TABLE IF EXISTS vlite_unique;
            DROP TABLE IF EXISTS catalog_match;
            DROP TABLE IF EXISTS corrected_flux;
            DROP TABLE IF EXISTS detected_source;
            DROP TABLE IF EXISTS detected_island;
            DROP TABLE IF EXISTS image;
            DROP TABLE IF EXISTS assoc_source;
            DROP TABLE IF EXISTS error;
            DROP TABLE IF EXISTS run_config;
            DROP FUNCTION IF EXISTS update_assoc_func;
            DROP FUNCTION IF EXISTS remove_vu_func;
            DROP FUNCTION IF EXISTS update_detected_func;
            ''')
        cur.execute(sql)

        print('\nCreating new tables...')
        sql = (
            '''
            CREATE EXTENSION IF NOT EXISTS q3c;

            CREATE TABLE run_config (
                id SERIAL NOT NULL,
                file TEXT,
                start_time TIMESTAMP (0),
                execution_time TIME (1),
                nimages INTEGER,
                stages JSON,
                options JSON,
                setup JSON,
                pybdsf_params JSON,
                image_qa_params JSON,
                PRIMARY KEY (id)
            );

            CREATE TABLE error (
                id INTEGER NOT NULL,
                reason TEXT,
                PRIMARY KEY (id)
            );

            CREATE TABLE assoc_source (
                id SERIAL NOT NULL,
                ra DOUBLE PRECISION,
                e_ra DOUBLE PRECISION,
                dec DOUBLE PRECISION,
                e_dec DOUBLE PRECISION,
                beam DOUBLE PRECISION,
                ndetect INTEGER,
                nmatches INTEGER,
                PRIMARY KEY (id)
            )
            WITH (fillfactor=90);

            CREATE TABLE image (
                id SERIAL NOT NULL UNIQUE,
                filename TEXT UNIQUE,
                imsize VARCHAR(14),
                obs_ra DOUBLE PRECISION,
                obs_dec DOUBLE PRECISION,
                pixel_scale DOUBLE PRECISION,
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
                mjdtime DOUBLE PRECISION,
                tau_time REAL,
                duration REAL,
                radius REAL,
                nsrc INTEGER,
                rms_box VARCHAR(14),
                stage INTEGER,
                catalogs_checked JSON,
                error_id INTEGER,
                nearest_problem TEXT,
                separation REAL,
                PRIMARY KEY (id),
                FOREIGN KEY (error_id) 
                  REFERENCES error (id) 
                  ON UPDATE CASCADE
            );

            CREATE TABLE detected_island (
                isl_id INTEGER NOT NULL,
                image_id INTEGER NOT NULL,
                total_flux DOUBLE PRECISION,
                e_total_flux DOUBLE PRECISION,
                rms DOUBLE PRECISION,
                mean DOUBLE PRECISION,
                resid_rms DOUBLE PRECISION,
                resid_mean DOUBLE PRECISION,
                PRIMARY KEY (isl_id, image_id),
                FOREIGN KEY (image_id) 
                    REFERENCES image (id) 
                    ON DELETE CASCADE
            );

            CREATE TABLE detected_source (
                src_id INTEGER NOT NULL,
                isl_id INTEGER NOT NULL,
                image_id INTEGER NOT NULL,
                ra DOUBLE PRECISION,
                e_ra DOUBLE PRECISION,
                dec DOUBLE PRECISION,
                e_dec DOUBLE PRECISION,
                total_flux DOUBLE PRECISION,
                e_total_flux DOUBLE PRECISION,
                peak_flux DOUBLE PRECISION,
                e_peak_flux DOUBLE PRECISION,
                ra_max DOUBLE PRECISION,
                e_ra_max DOUBLE PRECISION,
                dec_max DOUBLE PRECISION,
                e_dec_max DOUBLE PRECISION,
                maj DOUBLE PRECISION,
                e_maj DOUBLE PRECISION,
                min DOUBLE PRECISION,
                e_min DOUBLE PRECISION,
                pa DOUBLE PRECISION,
                e_pa DOUBLE PRECISION,
                dc_maj DOUBLE PRECISION,
                e_dc_maj DOUBLE PRECISION,
                dc_min DOUBLE PRECISION,
                e_dc_min DOUBLE PRECISION,
                dc_pa DOUBLE PRECISION,
                e_dc_pa DOUBLE PRECISION,
                code TEXT,
                assoc_id INTEGER,
                PRIMARY KEY (src_id, image_id),
                FOREIGN KEY (isl_id, image_id)
                  REFERENCES detected_island (isl_id, image_id)
                  ON DELETE CASCADE
            );

            CREATE TABLE corrected_flux (
                src_id INTEGER NOT NULL,
                isl_id INTEGER NOT NULL,
                image_id INTEGER NOT NULL,
                total_flux DOUBLE PRECISION,
                e_total_flux DOUBLE PRECISION,
                peak_flux DOUBLE PRECISION,
                e_peak_flux DOUBLE PRECISION,
                isl_total_flux DOUBLE PRECISION,
                isl_e_total_flux DOUBLE PRECISION,
                isl_rms DOUBLE PRECISION,
                isl_mean DOUBLE PRECISION,
                isl_resid_rms DOUBLE PRECISION,
                isl_resid_mean DOUBLE PRECISION,
                distance_from_center DOUBLE PRECISION,
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
                min_deRuiter DOUBLE PRECISION,
                PRIMARY KEY (id),
                FOREIGN KEY (assoc_id)
                  REFERENCES assoc_source (id)
                  ON DELETE CASCADE
            );

            CREATE TABLE vlite_unique (
                id SERIAL NOT NULL,
                image_id INTEGER,
                assoc_id INTEGER,
                detected BOOLEAN,
                PRIMARY KEY (id),
                FOREIGN KEY (image_id)
                  REFERENCES image (id)
                  ON DELETE CASCADE,
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

        make_error(cur, params)
        
        conn.commit()

        # Triggers
        sql = (
            '''
            CREATE OR REPLACE FUNCTION update_assoc_func() 
              RETURNS trigger AS $$
            BEGIN
              DELETE FROM assoc_source
              WHERE id = OLD.assoc_id AND ndetect = 1;
              UPDATE assoc_source SET 
                ra = (1./((1./(e_ra*e_ra))-(1./(OLD.e_ra*OLD.e_ra))))*(
                  (ra/(e_ra*e_ra))-(OLD.ra/(OLD.e_ra*OLD.e_ra))),
                e_ra = SQRT(1./((1./(e_ra*e_ra))-(1./(OLD.e_ra*OLD.e_ra)))),
                dec = (1./((1./(e_dec*e_dec))-(1./(OLD.e_dec*OLD.e_dec))))*(
                  (dec/(e_dec*e_dec))-(OLD.dec/(OLD.e_dec*OLD.e_dec))),
                e_dec = SQRT(1./(
                  (1./(e_dec*e_dec))-(1./(OLD.e_dec*OLD.e_dec)))),
                ndetect = ndetect - 1
              WHERE id = OLD.assoc_id;
            RETURN NEW;
            END;
            $$ LANGUAGE plpgsql;

            CREATE TRIGGER update_assoc
              AFTER DELETE ON detected_source
              FOR EACH ROW
              EXECUTE PROCEDURE update_assoc_func();
            ''')

        cur.execute(sql)

        sql = (
            '''
            CREATE OR REPLACE FUNCTION remove_vu_func()
              RETURNS TRIGGER AS $$
            BEGIN
              DELETE FROM vlite_unique WHERE assoc_id = OLD.id;
            RETURN NEW;
            END;
            $$ LANGUAGE plpgsql;

            CREATE TRIGGER remove_vu
              AFTER UPDATE OF nmatches ON assoc_source
              FOR EACH ROW
              WHEN (OLD.nmatches = 0 AND NEW.nmatches = 1)
              EXECUTE PROCEDURE remove_vu_func();
            ''')
        
        cur.execute(sql)

        sql = (
            '''
            CREATE OR REPLACE FUNCTION update_detected_func()
              RETURNS TRIGGER AS $$
            BEGIN
              UPDATE detected_source SET assoc_id = -1
              WHERE assoc_id = OLD.id;
            RETURN NEW;
            END;
            $$ LANGUAGE plpgsql;

            CREATE TRIGGER update_detected
              AFTER DELETE ON assoc_source
              FOR EACH ROW
              EXECUTE PROCEDURE update_detected_func();
            ''')

        cur.execute(sql)

        sql = (
            '''
            CREATE OR REPLACE FUNCTION update_nmatches_func()
              RETURNS TRIGGER AS $$
            DECLARE
              asid INTEGER;
              nm INTEGER;
            BEGIN
              UPDATE assoc_source SET nmatches = nmatches - 1
              WHERE id = OLD.assoc_id;
              SELECT INTO asid, nm id, nmatches FROM assoc_source
              WHERE id = OLD.assoc_id;
              IF nm = 0 THEN
                CREATE TABLE new_vu (
                  assoc_id INTEGER,
                  nmatches INTEGER
                );
                INSERT INTO new_vu (assoc_id, nmatches)
                VALUES (asid, nm);
              END IF;
            RETURN NEW;
            END;
            $$ LANGUAGE plpgsql;

            CREATE TRIGGER update_nmatches
              AFTER DELETE ON catalog_match
              FOR EACH ROW
              EXECUTE PROCEDURE update_nmatches_func();
            ''')
        
        cur.execute(sql)        
        conn.commit()
        cur.close()

    else:
        print('\nAborting... database left unchanged.')
        sys.exit(0)

