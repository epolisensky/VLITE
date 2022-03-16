"""Creates database tables, functions, and triggers."""
import sys
import logging


# create logger
createdb_logger = logging.getLogger('vdp.database.createdb')


def make_error(cur, params):
    """Inserts values into the database error table."""
    reason_dict = {'image manually unassociated, error unknown' : -1,
                   'image missing necessary header keyword(s)' : 1,
                   'number of visibilities < {}'.
                   format(params['min nvis']) : 2,
                   'sensitivity metric (noise x sqrt(int. time)) <= 0 or > {}'.
                   format(params['max sensitivity metric']): 3,
                   'beam axis ratio > {}'.
                   format(params['max beam axis ratio']) : 4,
                   'bad imaging target (NCP or planet)' : 5,
                   'problem source in image field-of-view' : 6,
                   'PyBDSF failed to process' : 7,
                   'zero sources extracted' : 8,
                   'source count metric > {}'.
                   format(params['max source count metric']) : 9,
                   'number of CLEAN iterations < {}'.
                   format(params['min niter']) : 10,
                   'bmin < {} or > {} pixels'.
                   format(params['min bpix'],params['max bpix']) : 11,
                   'image missing primary calibrators' : 12,
                   'image missing CLEAN components' : 13,
                   'image missing NX table' : 14}

    sql = 'INSERT INTO error (id, reason) VALUES (%s, %s);'
    for key, value in sorted(reason_dict.items(), key=lambda x: x[1]):
        cur.execute(sql, (value, key))


def create(conn, params, safe=False):
    """Creates new tables and triggers for the connected PostgreSQL 
    database by dropping tables if they exist. The current user
    must own the tables or be a superuser in order to drop them.
    USE WITH CAUTION! DELETING THE DATA CANNOT BE UNDONE.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    params : dict
        Dictionary of quality check requirements from the run
        configuration file.
    safe : bool, optional
        If ``False``, the user will be warned that existing data
        is about to be deleted and prompted to continue. Default
        value is ``False``.
    """
    if not safe:
        cont = input(('\nWARNING: Any existing tables and data in this '
                          'database will be deleted. Are you sure you want '
                          'to continue? '))
    else:
        cont = 'yes'
    if cont == 'y' or cont == 'yes':
        createdb_logger.info('Dropping tables if they exist...')
        cur = conn.cursor()
        sql = (
            '''
            DROP TABLE IF EXISTS new_vu;
            DROP TABLE IF EXISTS vlite_unique;
            DROP TABLE IF EXISTS catalog_match;
            DROP TABLE IF EXISTS corrected_flux;
            DROP TABLE IF EXISTS detected_source;
            DROP TABLE IF EXISTS detected_null;
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

        createdb_logger.info('Creating new tables...')
        sql = (
            '''
            CREATE EXTENSION IF NOT EXISTS q3c;
            CREATE EXTENSION IF NOT EXISTS quantile;

            CREATE TABLE run_config (
                id SERIAL NOT NULL,
                config_file TEXT,
                log_file TEXT,
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
                res_class VARCHAR(4),
                ndetect INTEGER,
                ns INTEGER,
                nc INTEGER,
                nm INTEGER,
                nmatches INTEGER,
                ave_total DOUBLE PRECISION,
                e_ave_total DOUBLE PRECISION,
                ave_peak DOUBLE PRECISION,
                e_ave_peak DOUBLE PRECISION,
                v_total DOUBLE PRECISION,
                v_peak DOUBLE PRECISION,
                eta_total DOUBLE PRECISION,
                eta_peak DOUBLE PRECISION,  
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
                config VARCHAR(3),
                cycle VARCHAR(3),
                semester VARCHAR(6),
                nvis INTEGER,
                niter INTEGER,
                mjdtime DOUBLE PRECISION,
                tau_time REAL,
                duration REAL,
                radius REAL,
                nsrc INTEGER,
                nclean INTEGER,
                rms_box VARCHAR(14),
                stage INTEGER,
                catalogs_checked JSON,
                error_id INTEGER,
                nearest_problem TEXT,
                separation REAL,
                glon DOUBLE PRECISION,
                glat DOUBLE PRECISION,
                az_star REAL,
                el_star REAL,
                pa_star REAL,
                az_end REAL,
                el_end REAL,
                pa_end REAL,
                az_i REAL,
                alt_i REAL,
                parang_i REAL,
                az_f REAL,
                alt_f REAL,
                parang_f REAL,
                lst_i REAL,
                lst_f REAL,
                pri_cals JSON,
                ass_flag BOOLEAN,
		nsn INTEGER,
                tsky REAL,
                square TEXT,
                sunsep REAL,
                pbkey TEXT,
                pb_flag BOOLEAN,
                ninterval INTEGER,
                max_dt REAL,
                nvisnx INTEGER,
                nbeam INTEGER,
                pbparangs REAL[],
                pbweights REAL[],
                pbza REAL[],
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

            CREATE TABLE detected_null (
                assoc_id INTEGER NOT NULL,
                image_id INTEGER NOT NULL,
                ra DOUBLE PRECISION,
                dec DOUBLE PRECISION,
                total_flux DOUBLE PRECISION,
                e_total_flux DOUBLE PRECISION,
                peak_flux DOUBLE PRECISION,
                e_peak_flux DOUBLE PRECISION,
                distance_from_center REAL,
                polar_angle REAL,
                snr REAL,
                PRIMARY KEY (assoc_id, image_id),
                FOREIGN KEY (image_id)
                  REFERENCES image (id)
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
                distance_from_center REAL,
                polar_angle REAL,
                snr REAL,
                compactness REAL,
                clean BOOLEAN,
                xpix REAL,
                ypix REAL,
                assoc_id INTEGER,
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
                separation REAL,
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
              WHEN (OLD.nmatches = 0 AND NEW.nmatches > 0)
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
                CREATE TABLE IF NOT EXISTS new_vu (
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
        createdb_logger.info('Aborting... database left unchanged.')
        sys.exit(0)

