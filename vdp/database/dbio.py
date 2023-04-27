"""Functions for specific I/O to the PostgreSQL database tables."""

import psycopg2
import psycopg2.extras
import json
import logging
from database import dbclasses
from math import sqrt
import numpy as np

# create logger
dbio_logger = logging.getLogger('vdp.database.dbio')


def weighted_avg_and_std(values, weights):
    """Return the weighted average and standard deviation.

    Parameters
    ----------
    values : numpy ndarray to be averaged
    weights : ndarray with value weights.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, sqrt(variance))



def record_config(conn, cfgfile, logfile, start_time, exec_time, nimages,
                  stages, opts, setup, sfparams, qaparams):
    """Inputs information about the current run of the pipeline
    into the database **run_config** table. All contents of the
    configuration file are stored as dictionaries/json.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    cfgfile : str
        Name of the YAML run configuration file.
    logfile : str
        Name of the log file.
    start_time : ``datetime.datetime`` instance
        Date & time at the start of the pipeline run.
    exec_time : ``datetime.timedelta`` instance
        Execution time of the pipeline run.
    nimages : int
        Number of images processed. (Technically,
        number of Image objects initialized.)
    stages : dict
        Dictionary of the *stages* section of the configuration file.
    opts : dict
        Dictionary of the *options* section of the configuration file.
    setup : dict
        Dictionary of the *setup* section of the configuration file.
    sfparams : dict
        Dictionary of the *pybdsf_params* section of the
        configuration file.
    qaparams : dict
        Dictionary of the *image_qa_params* section of the
        configuration file.
    """
    jstages = json.dumps(stages)
    jopts = json.dumps(opts)
    jsetup = json.dumps(setup)
    jsfparams = json.dumps(sfparams)
    jqaparams = json.dumps(qaparams)
    cur = conn.cursor()
    cur.execute('''INSERT INTO run_config (
        config_file, log_file, start_time, execution_time, nimages, 
        stages, options, setup, pybdsf_params, image_qa_params) VALUES (
        %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)''',
                (cfgfile, logfile, start_time, exec_time, nimages,
                 jstages, jopts, jsetup, jsfparams, jqaparams))
    conn.commit()
    cur.close()

       
def status_check(conn, impath):
    """Returns the id, highest completed stage, and radius used
    for source finding from the database **image** table or ``None``
    if the image filename is not in the database.

    """
    cur = conn.cursor()
    cur.execute('SELECT id, stage, radius FROM Image WHERE filename = %s',
                (impath, ))
    status = cur.fetchone()
    conn.commit()
    cur.close()

    return status


def add_image(conn, img, status, delete=False):
    """Inserts or updates rows in the database
    **image** table.
    
    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with attributes
        set from image header.
    status : tuple
        If ``None``, the image is new and is added to the
        database **image** table. If not ``None``, the image
        id stored in the tuple is used to update the correct
        row in the **image** table.
    delete : bool, optional
        If ``True``, rows in the **detected_island** database table
        with the appropriate 'image_id' will be deleted,
        cascading to the **detected_source** and **corrected_flux**
        tables and triggering updates on the **assoc_source**, 
        **catalog_match**, and **vlite_unique** tables. Default value
        is ``False``.

    Returns
    -------
    img : ``database.dbclasses.Image`` instance
        Initialized Image object with updated id attribute,
        if newly added to the database.
    """
    # Increase the image count by one
    dbclasses.Image.num_images += 1

    cur = conn.cursor()

    if img.pri_cals == None: jpri_cals = None
    else: jpri_cals=json.dumps(img.pri_cals)
    
    # Add new image to DB
    if status is None:
        dbio_logger.info('Adding new entry to image table.')
        sql = '''INSERT INTO image (
            filename, imsize, obs_ra, obs_dec, pixel_scale, object, obs_date, 
            map_date, obs_freq, primary_freq, bmaj, bmin, bpa, noise, peak, 
            config, cycle, semester, nvis, niter, mjdtime, tau_time, duration, 
            radius, nsrc, nclean, rms_box, stage, error_id, nearest_problem, 
            separation, glon, glat, lst_i, lst_f, az_star, el_star, pa_star, 
            az_end, el_end, pa_end, az_i, az_f, alt_i, alt_f, parang_i, 
            parang_f, pri_cals, ass_flag, nsn, tsky, square, sunsep, pbkey, 
            pb_flag, max_dt, nvisnx, ninterval, nbeam, pbparangs, pbweights, 
            pbza)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id'''
        vals = (img.filename, img.imsize, img.obs_ra, img.obs_dec, 
                img.pixel_scale, img.obj, img.obs_date, img.map_date,
                img.obs_freq, img.pri_freq, img.bmaj, img.bmin, img.bpa,
                img.noise, img.peak, img.config, img.cycle, img.semester,
                img.nvis, img.niter, img.mjdtime,
                img.tau_time, img.duration, img.radius, img.nsrc, img.nclean,
                img.rms_box, img.stage, img.error_id, img.nearest_problem,
                img.separation, img.glon, img.glat, img.lst_i, img.lst_f,
                img.az_star, img.el_star, img.pa_star, img.az_end, img.el_end,
                img.pa_end, img.az_i, img.az_f, img.alt_i, img.alt_f,
                img.parang_i, img.parang_f, jpri_cals, img.ass_flag, img.nsn,
                img.tsky, img.square, img.sunsep, img.pbkey, img.pb_flag,
                img.max_dt, img.nvisnx, img.ninterval, img.nbeam,
                np.array(img.pbparangs).tolist(),
                np.array(img.pbweights).tolist(),
                np.array(img.pbza).tolist())
        cur.execute(sql, vals)
        img.id = cur.fetchone()[0]
    # Update existing image entry
    else:
        dbio_logger.info('Updating existing entries in image table.')
        img.id = status[0]
        sql = '''UPDATE image SET filename = %s, imsize = %s, obs_ra = %s,
            obs_dec = %s, pixel_scale = %s, object = %s, obs_date = %s, 
            map_date = %s, obs_freq = %s, primary_freq = %s, bmaj = %s, 
            bmin = %s, bpa = %s, noise = %s, peak = %s, config = %s, 
            cycle = %s, semester = %s, nvis = %s, niter = %s, mjdtime = %s, 
            tau_time = %s, duration = %s, radius = %s, nsrc = %s, nclean = %s, 
            rms_box = %s, stage = %s, catalogs_checked = %s, error_id = %s, 
            nearest_problem = %s, separation = %s, glon = %s, glat = %s, 
            lst_i = %s, lst_f = %s, az_star = %s, el_star = %s, pa_star = %s, 
            az_end = %s, el_end = %s, pa_end = %s, az_i = %s, az_f = %s, 
            alt_i = %s, alt_f = %s, parang_i = %s, parang_f = %s, 
            pri_cals = %s, ass_flag = %s, nsn = %s, tsky = %s, square = %s, 
            sunsep = %s, pbkey = %s, pb_flag = %s, ninterval = %s, 
            max_dt = %s, nvisnx = %s, nbeam = %s, pbparangs = %s, 
            pbweights = %s, pbza = %s WHERE id = %s'''
        vals = (img.filename, img.imsize, img.obs_ra, img.obs_dec,
                img.pixel_scale, img.obj, img.obs_date, img.map_date,
                img.obs_freq, img.pri_freq, img.bmaj, img.bmin, img.bpa,
                img.noise, img.peak, img.config, img.cycle, img.semester,
                img.nvis, img.niter, img.mjdtime, img.tau_time, img.duration,
                img.radius, img.nsrc, img.nclean, img.rms_box, img.stage, None,
                img.error_id, img.nearest_problem, img.separation, img.glon,
                img.glat, img.lst_i, img.lst_f, img.az_star, img.el_star,
                img.pa_star, img.az_end, img.el_end, img.pa_end, img.az_i,
                img.az_f, img.alt_i, img.alt_f, img.parang_i, img.parang_f,
                jpri_cals, img.ass_flag, img.nsn, img.tsky, img.square,
                img.sunsep, img.pbkey, img.pb_flag, img.ninterval, img.max_dt,
                img.nvisnx, img.nbeam, np.array(img.pbparangs).tolist(),
                np.array(img.pbweights).tolist(),
                np.array(img.pbza).tolist(),img.id)
        cur.execute(sql, vals)
        if delete:
            # Delete corresponding sources
            dbio_logger.info('Removing previous sources...')
            cur.execute('DELETE FROM detected_island WHERE image_id = %s', (
                img.id, ))
            cur.execute('DELETE FROM vlite_unique WHERE image_id = %s', (
                img.id, ))

    conn.commit()
    cur.close()

    return img


def add_sources(conn, img, sources):
    """Inserts source finding and measurement results from
    PyBDSF stored as DetectedSource object attributes into
    the database **detected_island** and **detected_source**
    tables. The **image** table is also updated with some
    results from the source finding.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.    
    img : ``database.dbclasses.Image`` instance
        Image object with attributes updated with source
        finding results.
    sources : list
        DetectedSource objects whose attributes are the
        elliptical Gaussian fitting results.
    """
    cur = conn.cursor()

    # Update image table
    sql = '''UPDATE image SET radius = %s, nsrc = %s, nclean = %s, 
        rms_box = %s, error_id = %s, stage = %s WHERE id = %s'''
    vals = (img.radius, img.nsrc, img.nclean, img.rms_box,
            img.error_id, img.stage, img.id)
    cur.execute(sql, vals)

    if sources is not None:
        pass
    else:
        conn.commit()
        cur.close()
        return
    
    # Add sources to  detected_source and detected_island tables
    dbio_logger.info('Adding detected sources to database.')
    for src in sources:
        src.image_id = img.id
        # Insert values into detected_island table
        sql = '''INSERT INTO detected_island (
            isl_id, image_id, total_flux, e_total_flux, 
            rms, mean, resid_rms, resid_mean) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (isl_id, image_id)
            DO NOTHING'''
        vals = (src.isl_id, src.image_id, src.total_flux_isl,
                src.total_flux_islE, src.rms_isl, src.mean_isl,
                src.resid_rms, src.resid_mean)
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
        vals = (src.src_id, src.isl_id, src.image_id, src.ra, src.e_ra,
                src.dec, src.e_dec, src.total_flux,
                src.e_total_flux, src.peak_flux, src.e_peak_flux,
                src.ra_max, src.e_ra_max, src.dec_max,
                src.e_dec_max, src.maj, src.e_maj, src.min,
                src.e_min, src.pa, src.e_pa, src.dc_maj,
                src.e_dc_maj, src.dc_min, src.e_dc_min, src.dc_pa,
                src.e_dc_pa, src.code, src.assoc_id)
        cur.execute(sql, vals)

    conn.commit()
    cur.close()


def add_nulls(conn, sources):
    """Inserts null finding and measurement results from
    PyBDSF force-fitting stored as DetectedSource object 
    attributes into the database **detected_null** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.    
    sources : list
        DetectedSource objects whose attributes are the
        elliptical Gaussian force-fitting results.
    """
    cur = conn.cursor()

    if sources is not None:
        pass
    else:
        conn.commit()
        cur.close()
        return
    
    # Add sources to  detected_null table
    dbio_logger.info('Adding detected nulls to database.')
    for src in sources:
        # Insert values into detected_null table
        sql = '''INSERT INTO detected_null (
            assoc_id, image_id, ra, dec, total_flux, e_total_flux, 
            peak_flux, e_peak_flux, 
            distance_from_center, polar_angle, snr) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s)'''
        vals = (src.assoc_id, src.image_id, src.ra, src.dec, src.total_flux,
                src.e_total_flux, src.peak_flux, src.e_peak_flux,
                src.dist_from_center, src.polar_angle, src.snr)
        cur.execute(sql, vals)

    conn.commit()
    cur.close()
    


def add_corrected(conn, src, status=None):
    """Inserts primary beam corrected flux values into
    the database **corrected_flux** table.

    status : tuple
        If ``None``, a new src is added to the
        corrected_flux table. If not ``None``, the 
        src_id and image_id stored in the tuple is used to 
        update the correct row in the corrected_flux table.
    """
    cur = conn.cursor()

    #Add new:
    if status is None:
        cur.execute('''INSERT INTO corrected_flux (
            src_id, isl_id, image_id, total_flux, e_total_flux, peak_flux,
            e_peak_flux, isl_total_flux, isl_e_total_flux, isl_rms, isl_mean,
            isl_resid_rms, isl_resid_mean, distance_from_center, polar_angle, 
            snr, compactness, clean, xpix, ypix, assoc_id, nn_src_id, nn_dist) VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)''',
                    (src.src_id, src.isl_id, src.image_id, src.total_flux,
                     src.e_total_flux, src.peak_flux, src.e_peak_flux,
                     src.total_flux_isl, src.total_flux_islE, src.rms_isl,
                     src.mean_isl, src.resid_rms, src.resid_mean,
                     src.dist_from_center, src.polar_angle, src.snr,
                     src.compactness, src.clean, src.xpix, src.ypix,
                     src.assoc_id, src.nn_src_id, src.nn_dist))
    #Or update existing:
    else:
        src_id = status[0]
        img_id = status[1]
        sql = '''UPDATE corrected_flux SET total_flux = %s, e_total_flux = %s, peak_flux = %s,
            e_peak_flux = %s, isl_total_flux = %s, isl_e_total_flux = %s, isl_rms = %s, isl_mean = %s,
            isl_resid_rms = %s, isl_resid_mean = %s WHERE src_id = %s AND image_id = %s'''
        vals= (src.total_flux, src.e_total_flux, src.peak_flux, src.e_peak_flux,
                     src.total_flux_isl, src.total_flux_islE, src.rms_isl,
                     src.mean_isl, src.resid_rms, src.resid_mean, src_id, img_id)
        cur.execute(sql, vals)
        
    conn.commit()
    cur.close()


def add_assoc(conn, sources):
    """Adds a newly detected VLITE source to the
    **assoc_source** table and updates the 'assoc_id'
    for that source in the **detected_source** and
    **corrected_flux** tables.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        Newly detected VLITE sources stored as 
        DetectedSource objects.

    Returns
    -------
    sources : list
        DetectedSource objects with updated assoc_id
        attribute.
    """
    cur = conn.cursor()
    
    for src in sources:
        cur.execute('''INSERT INTO assoc_source (
            ra, e_ra, dec, e_dec, res_class, ndetect, ave_total, e_ave_total, 
            ave_peak, e_ave_peak, ns, nc, nm) VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id''',
                    (src.ra, src.e_ra, src.dec, src.e_dec,
                     src.res_class, src.ndetect, src.total_flux,
                     src.e_total_flux, src.peak_flux, src.e_peak_flux,
                     src.ns, src.nc, src.nm))
        src.id = cur.fetchone()[0]
        src.assoc_id = src.id
        cur.execute('''UPDATE detected_source SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))
        cur.execute('''UPDATE corrected_flux SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()
    
    return sources

def update_lcmetrics(conn, sources, flag=False):
    """Computes and updates the light curve metrics
    in the **assoc_source** table

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        DetectedSource objects extracted from the
        **assoc_source** table which have been matched
        to sources detected in the current image.
    flag : boolean
        If True, updates ave_peak, ave_total, e_ave_peak, e_ave_total
    """
    cur = conn.cursor()
    for src in sources:
        #fetch corrected fluxes and uncertainities
        cur.execute("SELECT SUM(peak_flux/(e_peak_flux*e_peak_flux)) AS peaktmp, SQRT(1./SUM(1./(e_peak_flux*e_peak_flux))) AS epeakflux, SUM(total_flux/(e_total_flux*e_total_flux)) AS totaltmp, SQRT(1./SUM(1./(e_total_flux*e_total_flux))) AS etotalflux FROM corrected_flux WHERE assoc_id = %s",(src.id,))
        arow=cur.fetchone()
        newpeak = arow[1]*arow[1]*arow[0]
        newtotal = arow[3]*arow[3]*arow[2]
        if flag:
            cur.execute("UPDATE assoc_source SET ave_peak = %s, e_ave_peak = %s, ave_total = %s, e_ave_total = %s WHERE id = %s",(newpeak, arow[1], newtotal, arow[3], src.id))
        #fetch light curve
        cur.execute('''SELECT total_flux, peak_flux, e_total_flux, e_peak_flux FROM corrected_flux WHERE assoc_id= %s''',(src.id,))
        numrows = int(cur.rowcount)
        #if only 1 detection then can't set metrics
        if numrows == 1:
            cur.execute("UPDATE assoc_source SET v_total = %s, v_peak = %s, eta_total = %s, eta_peak = %s WHERE id = %s",(None, None, None, None, src.id))
            #print('*',src.id,at,newtotal,ap,newpeak,' ',st,arow[3],sp,arow[1],' ',numrows)
        else:
            lc = np.fromiter(cur.fetchall(), dtype=[('total','float'),('peak','float'),
                            ('e_total','float'),('e_peak','float')], count=numrows)
            wt = 1./lc['e_total']**2
            wp = 1./lc['e_peak']**2
            # calc weighted average and std dev
            at,st = weighted_avg_and_std(lc['total'], wt)
            ap,sp = weighted_avg_and_std(lc['peak'], wp)
            # calc V
            src.v_total = st/at
            src.v_peak = sp/ap
            # calc eta
            src.eta_total = np.sum(wt*(lc['total']-at)**2)/(numrows-1)
            src.eta_peak  = np.sum(wp*(lc['peak']-ap)**2)/(numrows-1)
            # update
            cur.execute('''UPDATE assoc_source SET v_total = %s,
                v_peak = %s, eta_total = %s, eta_peak = %s WHERE id = %s''',
                (src.v_total, src.v_peak, src.eta_total, src.eta_peak, src.id))
            conn.commit()
            #print(' ',src.id,at,newtotal,ap,newpeak,' ',st,arow[3],sp,arow[1],' ',numrows)
    cur.close()

'''
def update_assoc(conn, sources):
    """Computes and updates the light curve metrics,
    positions and # detections in the **assoc_source** table

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        DetectedSource objects extracted from the
        **assoc_source** table which have been matched
        to sources detected in the current image.
    """
    cur = conn.cursor()

    for src in sources:
        #update counts
        cur.execute("SELECT COUNT(1) AS ndetect, SUM(CASE WHEN code='S' THEN 1 ELSE 0 END) AS ns, 
            SUM(CASE WHEN code='C' THEN 1 ELSE 0 END) AS nc, SUM(CASE WHEN code='M' THEN 1 ELSE 0 END) AS nm 
            FROM detected_source WHERE assoc_id = %s",(src.id,))
        row=cur.fetchone()
        cur.execute("UPDATE assoc_source SET ndetect = %s, ns = %s, nc = %s, nm = %s WHERE id = %s",
            (row['ndetect'], row['ns'], row['nc'], row['nm'], src.id))
        conn.commit()
        
        #fetch lightcurve: corrected fluxes and uncertainities
        cur.execute("SELECT total_flux, peak_flux, e_total_flux, e_peak_flux 
        FROM corrected_flux WHERE assoc_id= %s",(src.id,))
        numrows = int(cur.rowcount)
        lc = np.fromiter(cur.fetchall(), dtype=[('total','float'),('peak','float'),
                        ('e_total','float'),('e_peak','float')], count=numrows)
        wt = 1./lc['e_total']**2
        wp = 1./lc['e_peak']**2
        # calc weighted average and std dev
        at,st = weighted_avg_and_std(lc['total'], wt)
        ap,sp = weighted_avg_and_std(lc['peak'], wp)
        src.ave_peak = ap
        #src.e_ave_peak = sp
        src.ave_total = at
        #src.e_ave_total = st
        # calc V
        src.v_total = st/at
        src.v_peak = sp/ap
        # calc eta
        src.eta_total = np.sum(wt*(lc['total']-at)**2)/(numrows-1)
        src.eta_peak  = np.sum(wp*(lc['peak']-ap)**2)/(numrows-1)
        # update
        cur.execute("UPDATE assoc_source SET v_total = %s,
            v_peak = %s, eta_total = %s, eta_peak = %s WHERE id = %s",
            (src.v_total, src.v_peak, src.eta_total, src.eta_peak, src.id))

    conn.commit()
    cur.close()
'''
    
def update_matched_assoc(conn, sources):
    """Updates the RA & Dec of existing sources in the
    **assoc_source** table to the weighted average of all
    detections. Also updates the light curve average fluxes
    and number of detections

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        DetectedSource objects extracted from the
        **assoc_source** table which have been matched
        to sources detected in the current image.
    """
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE assoc_source SET ra = %s, e_ra = %s, dec = %s,
            e_dec = %s, ndetect = %s, ave_total = %s, e_ave_total = %s, 
            ave_peak = %s, e_ave_peak = %s, ns = %s, nc = %s, nm = %s 
            WHERE id = %s''',
                    (src.ra, src.e_ra, src.dec, src.e_dec, src.ndetect,
                     src.ave_total, src.e_ave_total, src.ave_peak,
                     src.e_ave_peak, src.ns, src.nc, src.nm, src.id))

    conn.commit()
    cur.close()


def update_detected_associd(conn, sources):
    """Updates the 'assoc_id' for sources in the
    **detected_source** table which have been successfully
    associated with existing VLITE sources in the 
    **assoc_source** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        DetectedSource objects detected in the current
        image that have been associated with previous
        VLITE sources in the **assoc_source** table.
    """
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE detected_source SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()

def update_corrected_associd(conn, sources):
    """Updates the 'assoc_id' for sources in the
    **corrected_flux** table which have been successfully
    associated with existing VLITE sources in the 
    **assoc_source** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        DetectedSource objects detected in the current
        image that have been associated with previous
        VLITE sources in the **assoc_source** table.
    """
    cur = conn.cursor()

    for src in sources:
        cur.execute('''UPDATE corrected_flux SET assoc_id = %s
            WHERE src_id = %s AND image_id = %s''',
                    (src.assoc_id, src.src_id, src.image_id))

    conn.commit()
    cur.close()


def update_checked_catalogs(conn, image_id, catalogs):
    """Updates the list of sky catalogs that have been
    checked for matches to VLITE sources detected in the
    image. A list of new catalogs to check is additionally
    defined by first querying the 'catalogs_checked' column
    in the **image** table and comparing to the list of catalogs
    specified from the configuration file.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    image_id : int
        Id number of the image in consideration.
    catalogs : list
        List of sky catalogs to use in matching as specified
        in the run configuration file.

    Returns
    -------
    new_catalogs : list
        Filtered list of only the sky catalogs for which
        cross-matching has not yet been carried out for
        the VLITE sources in the image.
    """
    cur = conn.cursor()

    # Get the catalogs which have already been checked for this image
    cur.execute('SELECT catalogs_checked FROM image WHERE id = %s',
                (image_id, ))
    already_checked = cur.fetchone()[0]
    if already_checked is not None:
        existing_catalogs = already_checked
    else:
        existing_catalogs = []

    # Filter out the already checked catalogs
    new_catalogs = [catalog for catalog in catalogs \
                    if catalog not in existing_catalogs]

    # Update catalogs_checked with new catalogs
    all_catalogs = json.dumps(sorted(existing_catalogs + new_catalogs))
    
    cur.execute('UPDATE image SET catalogs_checked = %s WHERE id = %s',
                (all_catalogs, image_id))

    conn.commit()
    cur.close()

    return new_catalogs


def check_catalog_match(conn, asrc_id, catalog):
    """Checks if a source from the **assoc_source** table 
    already has a match to a sky catalog source in the
    **catalog_match** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    asrc_id : int
        Row id of the source in the **assoc_soure** table to
        match to the 'assoc_id' in the **catalog_match** table.
    catalog : str
        Name of the sky catalog.

    Returns
    -------
    rowid : int
        The row id of the entry if it exists. Returns
        ``None`` otherwise.
    """
    cur = conn.cursor()

    cur.execute('''SELECT id FROM catalog_match
        WHERE assoc_id = %s AND catalog_id = (
          SELECT id FROM radcat.catalogs WHERE name = %s)''',
                (asrc_id, catalog))
    rowid = cur.fetchone()

    conn.commit()
    cur.close()

    return rowid


def update_assoc_nmatches(conn, sources):
    """Updates the number of sky catalog matches to
    a given source in the **assoc_source** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        DetectedSource objects with updated nmatches
        attribute after running sky catalog matching,
        or list of assoc_ids.
    """
    cur = conn.cursor()

    for src in sources:
        try:
            # sources = list of objects
            src_id = src.id
            nmatches = src.nmatches
        except AttributeError:
            # sources = list of assoc_ids
            src_id = src
            cur.execute('SELECT nmatches FROM assoc_source WHERE id = %s',
                        (src_id, ))
            nmatches = cur.fetchone()[0] + 1
        cur.execute('UPDATE assoc_source SET nmatches = %s WHERE id = %s',
                    (nmatches, src_id))

    conn.commit()
    cur.close()


def add_catalog_match(conn, sources):
    """Adds an entry to the **catalog_match** table for
    every sky catalog source matched to a VLITE source
    in the **assoc_source** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        List of CatalogSource objects matched to VLITE
        DetectedSource objects, or list of tuples.
    """
    cur = conn.cursor()

    for src in sources:
        try:
            # sources = list of objects
            catalog_id = src.catalog_id
            src_id = src.id
            assoc_id = src.assoc_id
            sep = src.sep
        except AttributeError:
            # sources = list of tuples
            cur.execute('SELECT id FROM radcat.catalogs WHERE name = %s',
                        (src[0], ))
            catalog_id = cur.fetchone()[0]
            src_id = src[1]
            assoc_id = src[2]
            sep = src[3]
        cur.execute('''INSERT INTO catalog_match (
            catalog_id, src_id, assoc_id, separation) VALUES (
            %s, %s, %s, %s)''',
                    (catalog_id, src_id, assoc_id, sep))

    conn.commit()
    cur.close()


def check_vlite_unique(conn, asrc_id):
    """Checks if a given source from the **assoc_source**
    table is already in the **vlite_unique** table. This
    is so that sources don't get added twice to the
    **vlite_unique** table (once when the nmatches = 0 source
    is pulled from **assoc_source** table and again if no
    sky catalog match is found) when updating the catalog
    matching results by adding new sky catalogs without
    re-doing the previous catalog matching results.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    asrc_id : int
        Row id of the source in the **assoc_soure** table to
        match to the 'assoc_id' in the **vlite_unique** table.

    Returns
    -------
    existing : list
        Returns the 'id', 'image_id', and 'detected' columns of
        the entry with the given 'assoc_id'. Otherwise, 
        returns ``None``.
    """
    cur = conn.cursor()

    cur.execute('''SELECT id, image_id, detected FROM vlite_unique
        WHERE assoc_id = %s''', (asrc_id, ))
    existing = cur.fetchall()

    conn.commit()
    cur.close()
   
    return existing


def add_vlite_unique(conn, src, image_id, update=False):
    """Adds or updates an entry in the **vlite_unique** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    src : ``database.dbclasses.DetectedSource`` instance 
        DetectedSource object with attribute nmatches = 0.
    image_id : int
        Id number of the image from which the source came.
    update : bool, optional
        If ``True``, the 'detected' column is updated for
        the existing row with the specified 'image_id' and
        'assoc_id'. Otherwise, a new row is added.
        Default is ``False``.
    """
    cur = conn.cursor()

    if update:
        cur.execute('''UPDATE vlite_unique SET detected = %s
            WHERE image_id = %s AND assoc_id = %s''',
                    (src.detected, image_id, src.id))
    else:
        cur.execute('''INSERT INTO vlite_unique (
            image_id, assoc_id, detected) VALUES (%s, %s, %s)''',
                    (image_id, src.id, src.detected))

    conn.commit()
    cur.close()


def get_image_sources(conn, image_id):
    """Returns a list of sources belonging to a
    particular image from the **detected_source** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    image_id : int
        Id number of the image.

    Returns
    -------
    detected_sources : list
        Sources pulled from the **detected_source** table
        translated from row dictionary objects to
        DetectedSource objects.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    cur.execute('SELECT * FROM detected_source WHERE image_id = %s',
                (image_id, ))
    rows = cur.fetchall()

    conn.commit()

    detected_sources = []
    for row in rows:
        detected_sources.append(dbclasses.DetectedSource())
        dbclasses.dict2attr(detected_sources[-1], row)
        cur.execute('SELECT * FROM detected_island WHERE isl_id = %s AND image_id = %s', (detected_sources[-1].isl_id,detected_sources[-1].image_id))
        row2=cur.fetchone()
        detected_sources[-1].rms_isl         = row2['rms']
        detected_sources[-1].mean_isl        = row2['mean']
        detected_sources[-1].resid_rms       = row2['resid_rms']
        detected_sources[-1].resid_mean      = row2['resid_mean']
        detected_sources[-1].total_flux_isl  = row2['total_flux']
        detected_sources[-1].total_flux_islE = row2['e_total_flux']

    cur.close()
    return detected_sources


def get_associated(conn, sources):
    """Returns a list of sources belonging to a
    particular image from the **assoc_source** table.
    
    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        TWO OPTIONS:
        1) DetectedSource objects pulled from the
        **detected_source** table based on 'image_id'.
        2) List of assoc_source ids 

    Returns
    -------
    assoc_sources : list
        Sources pulled from the **assoc_source** table
        based on matching 'assoc_id' and translated
        from row dictionary objects to 
        DetectedSource objects
    """   
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    rows = []
    for src in sources:
        try:
            # OPTION 1: sources = list of objects 
            src_id = src.assoc_id
        except AttributeError:
            # OPTION 2: sources = list of assoc_ids
            src_id = src
        cur.execute('SELECT * FROM assoc_source WHERE id = %s',
                    (src_id, ))
        rows.append(cur.fetchone())

    conn.commit()
    cur.close()

    assoc_sources = []
    for row in rows:
        if not row:
            continue
        assoc_sources.append(dbclasses.DetectedSource())
        dbclasses.dict2attr(assoc_sources[-1], row)

    return assoc_sources


def delete_matches(conn, sources, image_id):
    """Deletes all previous sky catalog cross-matching
    results for a given set of sources. This function
    is called when the configuration file option *redo match*
    is ``True``.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        DetectedSource objects belonging to a
        particular image pulled from the
        **assoc_source** table based on matching 'assoc_id'.
    image_id : int
        Id number of the image.

    Returns
    -------
    sources : list
        DetectedSource objects with their nmatches
        attribute re-initialized to 0.
    """
    dbio_logger.info('Removing previous sky catalog matching results '
                     'for {} sources.'.format(len(sources)))

    cur = conn.cursor()

    cur.execute('UPDATE image SET catalogs_checked = %s WHERE id = %s',
                (None, image_id))

    for src in sources:
        src.nmatches = 0
        cur.execute('UPDATE assoc_source SET nmatches = %s WHERE id = %s',
                    (src.nmatches, src.id))
        cur.execute('DELETE FROM catalog_match WHERE assoc_id = %s',
                    (src.id, ))
        cur.execute('''DELETE FROM vlite_unique WHERE image_id = %s
            AND assoc_id = %s''', (image_id, src.id))

    conn.commit()
    cur.close()

    return sources


def update_stage(conn, imobj):
    """Updates the stage column in the **image** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Image object used to set table column values.
    """
    cur = conn.cursor()

    cur.execute('UPDATE image SET stage = %s WHERE id = %s',
                (imobj.stage, imobj.id))

    conn.commit()
    cur.close()


def remove_catalog(conn, catalogs):
    """Deletes rows with results from specified catalogs
    from the **catalog_match** table.

    """
    cur = conn.cursor()

    for catalog in catalogs:
        # Get the catalog_id from the radcat.catalogs table
        cur.execute('SELECT id FROM radcat.catalogs WHERE name = %s',
                    (catalog, ))
        catid = cur.fetchone()
        if catid is None:
            dbio_logger.error('ERROR: the catalog {} does not exist.'.
                              format(catalog))
            return
        else: pass

        # Nothing will happen if the catalog_id is not in the table
        cur.execute('DELETE FROM catalog_match WHERE catalog_id = %s',
                    (catid[0], ))
        # Find all images which need their catalogs_checked list updated
        cur.execute('''SELECT id, catalogs_checked FROM image
            WHERE catalogs_checked::jsonb ? %s''',
                    (catalog, ))
        rows_to_update = cur.fetchall()
        
        updated_rows = []
        if rows_to_update is not None:
            for row in rows_to_update:
                old_catalogs = row[1]
                new_catalogs = [cat for cat in old_catalogs \
                                if cat not in [catalog]]
                updated_rows.append((row[0], json.dumps(new_catalogs)))

        if updated_rows:
            for row in updated_rows:
                cur.execute('''UPDATE image SET catalogs_checked = %s
                    WHERE id = %s''', (row[1], row[0]))

        conn.commit()
        
    cur.close()


def get_new_vu(conn):
    """Returns list of sources from the **assoc_source** table
    for which the 'nmatches' dropped to 0 after removing 
    catalog matching results. The **new_vu** table
    is created to record the row id number of sources in the
    **assoc_source** table whose 'nmatches' column is 0 after
    subtracting 1 as triggered by deletion of rows in the
    **catalog_match** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.

    Returns
    -------
    assoc_objs : list
        List of DetectedSource objects pulled from the
        **assoc_source** table whose 'nmatches' now equals 0.
    """
    cur = conn.cursor()

    # Get assoc_ids of all sources whose nmatches went down to 0
    try:
        cur.execute('SELECT assoc_id FROM new_vu')
        asrc_ids = cur.fetchall()
        cur.execute('DROP TABLE new_vu')
    except psycopg2.ProgrammingError:
        asrc_ids = None

    conn.commit()
    cur.close()

    # No new_vu table; no new VU sources
    if asrc_ids is None:
        return None

    # Create list of objects for all assoc_sources which are new VU sources
    assoc_objs = get_associated(conn, asrc_ids)

    return assoc_objs


def get_vu_image(conn, assoc_id):
    """Retrieves the id numbers and field-of-view radii for
    any image which contains the specified the VLITE source.
    This information is needed to update the **vlite_unique**
    table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    assoc_id : int
        The id number of the VLITE source from the
        **assoc_source** table.

    Returns
    -------
    vu_image_list : list
        List of tuples containing the id number and radius of
        any image in which the specified VLITE source was
        detected.
    """
    cur = conn.cursor()

    cur.execute('SELECT image_id FROM detected_source WHERE assoc_id = %s',
                (assoc_id, ))
    rows = cur.fetchall()
    image_ids = [row[0] for row in rows]
    
    radii = []
    for image_id in image_ids:
        cur.execute('SELECT radius FROM image WHERE id = %s',
                    (image_id, ))
        radii.append(cur.fetchone()[0])

    conn.commit()
    cur.close()

    vu_image_list = zip(image_ids, radii)

    return vu_image_list


def remove_sources(conn, assoc_ids):
    """Deletes the specified sources from the database
    **assoc_source** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    assoc_ids : list
        List of id numbers corresponding to the VLITE sources
        to be removed from the **assoc_source** table.    
    """
    cur = conn.cursor()

    cur.execute('DELETE FROM assoc_source WHERE id IN %s',
                (assoc_ids, ))

    conn.commit()
    cur.close()


def remove_images(conn, images):
    """Deletes the specified images from the database.
    Removal propagates to all affected tables.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    images : list
        List of image filenames to be removed from
        the database **image** table.
    """
    cur = conn.cursor()

    for image in images:
        cur.execute('SELECT id FROM image WHERE filename LIKE %s',('%'+image, ))
        try:
            image_id = cur.fetchone()[0]
        except TypeError:
            dbio_logger.info('WARNING: Image {} is not in the database.'.
                             format(image))
        #fetch assoc_id's for sources w/ ndetect >= 2 in this image
        cur.execute("SELECT a.id FROM assoc_source AS a, detected_source AS d WHERE a.id=d.assoc_id AND d.image_id = %s AND a.ndetect > 1 ORDER BY a.id",(image_id,))
        assid = cur.fetchall()
        #get assoc_source objects
        asrcs=get_associated(conn, assid)
        #delete image
        cur.execute('DELETE FROM image WHERE id = %s',(image_id, ))
        conn.commit()
        #update fluxes, LC metrics for assoc_sources
        update_lcmetrics(conn,asrcs,True)
    cur.close()


def unassoc_images(conn, images):
    """Unassociates the specified images from the database.
    Updates image, assoc_source, detected_source, corrected_flux tables.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    images : list
        List of image filenames to be unassociated from
        the database.
    """
    cur = conn.cursor()

    for image in images:
        #print('image = ',image)
        cur.execute('SELECT id,stage,filename FROM image WHERE filename LIKE %s',
                    ('%'+image, ))
        try:
            row = cur.fetchone()
            image_id = row[0]
            stage = int(row[1])
            #print(image_id,stage)
            if stage > 2:
                #fetch distinct assoc_id's for sources in this image
                cur.execute("SELECT DISTINCT(assoc_id) FROM detected_source WHERE image_id = %s ORDER BY assoc_id",(image_id,))
                assid = cur.fetchall()
                #update 
                cur.execute("UPDATE image SET stage = 2, error_id = %s WHERE id = %s",(-1,image_id,))
                cur.execute("UPDATE detected_source SET assoc_id = %s WHERE image_id = %s",(None,image_id))
                cur.execute("UPDATE corrected_flux SET assoc_id = %s WHERE image_id = %s",(None,image_id))
                conn.commit()
                #update values for each assoc_source in assid
                assoc_srcs = get_associated(conn,assid)
                for src in assoc_srcs:
                    #check if no detections left & delete
                    cur.execute("SELECT COUNT(1) FROM detected_source WHERE assoc_id = %s",(src.id,))
                    count=cur.fetchone()
                    if count[0]==0:
                        dbio_logger.info('Deleting assoc_source id {}'.format(src.id))
                        cur.execute('DELETE FROM assoc_source WHERE id = %s',(src.id,))
                        conn.commit()
                    else:
                        dbio_logger.info('Updating assoc_source id {}'.format(src.id))
                        #update coords, fluxes, # detections
                        cur.execute("SELECT SUM(d.ra/(d.e_ra*d.e_ra)) AS ratmp, SQRT(1./SUM(1./(d.e_ra*d.e_ra))) AS era, SUM(d.dec/(d.e_dec*d.e_dec)) AS dectmp, SQRT(1./SUM(1./(d.e_dec*d.e_dec))) AS edec, SUM(c.peak_flux/(c.e_peak_flux*c.e_peak_flux)) AS peaktmp, SQRT(1./SUM(1./(c.e_peak_flux*c.e_peak_flux))) AS epeakflux, SUM(c.total_flux/(c.e_total_flux*c.e_total_flux)) AS totaltmp, SQRT(1./SUM(1./(c.e_total_flux*c.e_total_flux))) AS etotalflux, COUNT(1) AS ndetect, SUM(CASE WHEN d.code='S' THEN 1 ELSE 0 END) AS ns, SUM(CASE WHEN d.code='C' THEN 1 ELSE 0 END) AS nc, SUM(CASE WHEN d.code='M' THEN 1 ELSE 0 END) AS nm FROM detected_source AS d, corrected_flux AS c WHERE (c.src_id,c.image_id,d.assoc_id) = (d.src_id,d.image_id,%s)",(src.id,))
                        arow=cur.fetchone()
                        newra = arow[1]*arow[1]*arow[0]
                        newdec = arow[3]*arow[3]*arow[2]
                        newpeak = arow[5]*arow[5]*arow[4]
                        newtotal = arow[7]*arow[7]*arow[6]
                        #print(newra,newdec,newpeak,newtotal,arow[8])
                        cur.execute("UPDATE assoc_source SET ra = %s, e_ra = %s, dec = %s, e_dec = %s, ave_peak = %s, e_ave_peak = %s, ave_total = %s, e_ave_total = %s, ndetect = %s, ns = %s, nc = %s, nm = %s WHERE id = %s",(newra, arow[1], newdec, arow[3], newpeak, arow[5], newtotal, arow[7], arow[8], arow[9], arow[10], arow[11], src.id))
                        conn.commit()
                        #update LC metrics
                        #print('updating LC metrics...')
                        #fetch lightcurve: corrected fluxes and uncertainities
                        cur.execute("SELECT total_flux, peak_flux, e_total_flux, e_peak_flux FROM corrected_flux WHERE assoc_id= %s",(src.id,))
                        numrows = int(cur.rowcount)
                        #print(numrows)
                        if numrows == 1:
                            #print('not enough pts for metrics')
                            cur.execute("UPDATE assoc_source SET v_total = %s, v_peak = %s, eta_total = %s, eta_peak = %s WHERE id = %s",
                                        (None, None, None, None, src.id))
                        else:
                            lc = np.fromiter(cur.fetchall(), dtype=[('total','float'),('peak','float'),
                                                                ('e_total','float'),('e_peak','float')], count=numrows)
                            wt = 1./lc['e_total']**2
                            wp = 1./lc['e_peak']**2
                            # calc weighted average and std dev
                            at,st = weighted_avg_and_std(lc['total'], wt)
                            ap,sp = weighted_avg_and_std(lc['peak'], wp)
                            # calc V
                            src.v_total = st/at
                            src.v_peak = sp/ap
                            # calc eta
                            src.eta_total = np.sum(wt*(lc['total']-at)**2)/(numrows-1)
                            src.eta_peak  = np.sum(wp*(lc['peak']-ap)**2)/(numrows-1)
                            # update
                            #print('  ',(src.v_total, src.v_peak, src.eta_total, src.eta_peak))
                            cur.execute("UPDATE assoc_source SET v_total = %s, v_peak = %s, eta_total = %s, eta_peak = %s WHERE id = %s",
                                        (src.v_total, src.v_peak, src.eta_total, src.eta_peak, src.id))
                        conn.commit()

        except TypeError:
            dbio_logger.info('WARNING: Image {} is not in the database.'.
                             format(image))
    cur.close()

    
def update_pbflag(conn, imobj):
    """Updates the pb_flag column in the **image** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Image object used to set table column values.
    """
    cur = conn.cursor()

    cur.execute('UPDATE image SET pb_flag = %s WHERE id = %s',
                (imobj.pb_flag, imobj.id))

    conn.commit()
    cur.close()

def update_assflag(conn, imobj):
    """Updates the ass_flag column in the **image** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Image object used to set table column values.
    """
    cur = conn.cursor()

    cur.execute('UPDATE image SET ass_flag = %s WHERE id = %s',
                (imobj.ass_flag, imobj.id))

    conn.commit()
    cur.close()

def update_pbvalues(conn, imobj):
    """Updates the pbparangs & pbweights column in the **image** table.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    imobj : ``database.dbclasses.Image`` instance
        Image object used to set table column values.
    """
    cur = conn.cursor()

    cur.execute('UPDATE image SET pbparangs = %s, pbweights = %s, pbza = %s WHERE id = %s',
                (np.array(imobj.pbparangs).tolist(), np.array(imobj.pbweights).tolist(), np.array(imobj.pbza).tolist(), imobj.id))
    if imobj.vcss:
        cur.execute('UPDATE image SET nbeam = %s WHERE id = %s',
                    (imobj.nbeam, imobj.id))

    conn.commit()
    cur.close()
    
''' NEEDS UPDATING
def update_corrected(conn,setup):
    """Updates the **corrected_flux** table with primary beam corrections.
    Reads image table, fetches all detected sources for each image,
    corrects source fluxes for primary beam.
    
    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute('SELECT id,filename,imsize,obs_ra,obs_dec,pixel_scale,object,obs_date,map_date,obs_freq,primary_freq,bmaj,bmin,bpa,noise,peak,config,nvis,niter,mjdtime,tau_time,duration,nsrc,rms_box,error_id,stage,radius,nearest_problem,separation FROM image')
    rows = cur.fetchall()
    cur.close()

    for row in rows:
        print(row)
        imobj=dbclasses.Image()
        dbclasses.dict2attr(imobj, row)
        imobj.pri_freq = row['primary_freq']
        img_id=imobj.id
        priband=imobj.pri_freq
        sources=get_image_sources(conn,img_id)
        for src in sources:
            src.calc_center_dist(imobj)
            src.correct_flux(priband)
            add_corrected(conn,src,status=(src.src_id,src.image_id))

    conn.commit()
'''
