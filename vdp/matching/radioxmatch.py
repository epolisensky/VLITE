"""radioxmatch.py contains all the machinery for cross-matching
radio point sources.

"""
import re
import numpy as np
import logging
import psycopg2
import psycopg2.extras
from psycopg2 import sql
from radiocatalogs import catalogio
from database import dbclasses, dbio
from matching import matchfuncs
from sklearn.neighbors import BallTree

# create logger
match_logger = logging.getLogger('vdp.matching.radioxmatch')


def write_regions(srclist, impath, ext='.reg'):
    """Writes a ds9 regions file from given source list.
    
    Parameters
    ----------
    srclist : list
        List of DetectedSource or CatalogSource objects
        to write to a regions file.
    impath : str
        Name of the image including full directory path.
    ext : str, optional
        Extension for the file name
        (i.e. '1.5GHz.0137+331.IPln1.matches.reg').
        Default is '.reg'.
    """
    fname = impath[:-5] + ext
    with open(fname, 'w') as f:
        f.write('global color=cyan font="helvetica 10 normal" '
                'select=1 highlite=1 edit=1 move=1 delete=1 '
                'include=1 fixed=0 source\n')
        f.write('fk5\n')
        for src in srclist:
            f.write('ellipse(%f,%f,%.2f",%.2f",%.1f) # text={%s}\n' % (
                src.ra, src.dec, src.maj, src.min, src.pa + 90.0, src.name))


def filter_res(rows, res_class):
    """Filters out sources extracted from the database
    which originate from images with spatial resolutions
    outside the acceptable range.

    Parameters
    ----------
    rows : list
        psycopg2 row dictionary objects extracted
        from the database.
    res_class : float
        Spatial resolution of the current image in arcseconds.

    Returns
    -------
    keep : list
        List of psycopg2 row dictionary objects after
        applying the spatial resolution filtering.
    """
    keep = []
    match_logger.info('Limiting to sources in resolution class {} '.format(res_class))

    for row in rows:
        if row['res_class'] == res_class:
            keep.append(row)

    match_logger.info(' -- {} sources remaining'.format(len(keep)))

    return keep


def check_previous(conn, src, search_radius):
    """Searches the database **image** table for images
    of the same size which could have contained the source
    in question in their field-of-view.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    src : ``database.dbclasses.DetectedSource`` instance
        Source whose position will be used to find
        all previously processed images which could
        have contained it.
    search_radius : float
        Radius of the image field-of-view in degrees.

    Returns
    -------
    prev_images : list
        List of ids of the images which could have
        contained the given point.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute('''SELECT id FROM image 
        WHERE q3c_join(%s, %s, obs_ra, obs_dec, %s) AND
        radius = %s''',
                (src.ra, src.dec, search_radius, search_radius))
    prev_images = cur.fetchall()
    cur.close()

    return prev_images


def cone_search(conn, table, center_ra, center_dec, radius, schema='public'):
    """Extracts all sources from the specified database
    table which fall within a circular region on the sky.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    table : str
        Name of the database table to query.
    center_ra : float
        Right ascension coordinate of the circle center
        in degrees.
    center_dec : float
        Declination coordinate of the circle center in degrees.
    radius : float
        Size of the circular search region in degrees.
    schema : str, optional
        The database schema which contains the table.
        Default is 'public'.

    Returns
    -------
    rows : list
        List of row dictionary objects corresponding to the
        sources that fall within the circular search region.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    query = 'SELECT * FROM {}.{} WHERE q3c_join(%s, %s, ra, dec, %s)'
    # This one isn't needed until > 1 million rows
    #query = 'SELECT * FROM {}.{} WHERE q3c_radial_query(ra, dec, %s, %s, %s)'
    cur.execute(psycopg2.sql.SQL(query).format(
        psycopg2.sql.Identifier(schema), psycopg2.sql.Identifier(table)),
                (center_ra, center_dec, radius))
    rows = cur.fetchall()
    cur.close()

    return rows


def associate(conn, detected_sources, imobj, search_radius, save):
    """Associates new sources with old sources if the center
    positions of the two sources are separated by an angular
    distance less than the full size of the major axis of
    the current image's beam.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    detected_sources : list
        List of the sources extracted from the current image.
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    search_radius : float
        Size of the circular search region in degrees.
    save : bool
        If ``False``, print the results to the console and/or
        log file.

    Returns
    -------
    detected_matched : list
        Sources extracted from the image that were successfully
        associated with previously detected sources stored in the
        **assoc_source** table.
    detected_unmatched : list
        Sources extracted from the image which could NOT be successfully
        associated with previously detected sources. These are added
        to the **assoc_table** as new detections and then cross-matched
        with other radio catalogs.
    assoc_matched : list
        Previously detected sources from the **assoc_source** table which
        were successfully associated with sources from the new image. 
        Positions of these sources are updated in the **assoc_source**
        table to reflect the weighted average of all detections. The 
        number of detections ('ndetect') is increased by one.
    assoc_unmatched : list
        Previously detected sources from the **assoc_source** table which
        were not successfully associated with sources from the new image.
        Any of these non-detections with no radio catalog matches
        ('nmatches' = 0) are recorded in the **vlite_unique** table.
    """
    if imobj.vcss or imobj.mosaic:
        res_class = 'VCSS'
    else:
        res_class = imobj.config

    # Extract all previously detected sources in the same FOV
    assoc_rows = cone_search(conn, 'assoc_source', imobj.obs_ra,
                             imobj.obs_dec, search_radius)
    match_logger.info('Extracted {} sources from assoc_source table '
                      'within {} degrees.'.format(
                          len(assoc_rows), search_radius))
    # Limit to sources taken from images of similar resolution
    if len(assoc_rows) > 0:
        filtered_assoc_rows = filter_res(assoc_rows, res_class)
    else:
        filtered_assoc_rows = []


    #for dsrc in detected_sources:
    dra=[]
    ddec=[]
    for i in range(len(detected_sources)):
        dra.append(detected_sources[i].ra)
        ddec.append(detected_sources[i].dec)
        #print(' ',detected_sources[i].src_id,detected_sources[i].ra,detected_sources[i].dec)
    #    print('  ',dsrc.src_id,dsrc.ra,dsrc.dec)
    #print(' ',detected_sources['src_id'],detected_sources['ra'],detected_sources['dec'])
    dra=np.array(dra)
    ddec=np.array(ddec)

    if not filtered_assoc_rows:
        # No previous sources found in that sky region at that resolution
        for src in detected_sources:
            src.res_class = res_class
            src.ndetect = 1
        detected_matched = []
        detected_unmatched = detected_sources
        assoc_matched = []
        assoc_unmatched = []
    else:
        # Translate row dictionaries to DetectedSource objects
        assoc_sources = []
        assoc_ids = []
        ara = []
        adec = []
        for asrc in filtered_assoc_rows:
            ara.append(asrc['ra'])
            adec.append(asrc['dec'])
            assoc_ids.append(asrc['id'])
            assoc_sources.append(dbclasses.DetectedSource())
            dbclasses.dict2attr(assoc_sources[-1], asrc)
        match_logger.info('Attempting to match {} sources from this image to '
                          '{} sources previously detected in VLITE images...'.
                          format(len(detected_sources), len(assoc_sources)))
        ara=np.array(ara)
        adec=np.array(adec)

        detected_matched = []
        detected_unmatched = []
        assoc_matched = []
        assoc_unmatched = []

        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        # Print results without saving to database
        if not save:
            # Dump detected_sources into temporary table
            sql = (
                '''
                CREATE TEMP TABLE temp_source (
                    src_id INTEGER,
                    ra DOUBLE PRECISION,
                    dec DOUBLE PRECISION
                );
                ''')
            cur.execute(sql)
            conn.commit()
            for src in detected_sources:
                cur.execute('''INSERT INTO temp_source (
                    src_id, ra, dec) VALUES (%s, %s, %s)''', (
                        src.src_id, src.ra, src.dec))
            conn.commit()
            # Find nearest neighbor & "match" if within half a beam
            sql = '''SELECT a.src_id, bb.id AS assoc_id,
                3600*q3c_dist(a.ra, a.dec, bb.ra, bb.dec) AS sep,
                3600*q3c_dist(a.ra, a.dec, bb.ra, bb.dec) < %s AS match
                FROM temp_source AS a, LATERAL (
                SELECT b.* FROM assoc_source AS b WHERE b.id IN %s
                ORDER BY q3c_dist(a.ra, a.dec, b.ra, b.dec) ASC LIMIT 1)
                AS bb'''
            values = (imobj.bmaj, tuple(assoc_ids))
            cur.execute(sql, values)
            rows = cur.fetchall()
            cur.execute('DROP TABLE temp_source')
            conn.commit()
            match_logger.info('-----------------------------------------------'
                              '-----------------------------------------------'
                              '---------------------------------')
            match_logger.info('src_id match assoc_id\tra\t\te_ra\t\t\tdec\t\t'
                              'e_dec\t\tseparation (arcsec)\tndetect')
            match_logger.info('-----------------------------------------------'
                              '-----------------------------------------------'
                              '---------------------------------')
        # Save association results for database
        else:
            # Find nearest neighbor & "match" if within half a beam
            # Make BallTree with assoc sources
            #tree = BallTree(np.column_stack((np.radians(assoc_sources['dec']), np.radians(assoc_sources['ra']))), leaf_size=20, metric='haversine')
            #returned is the nearest assoc neighbor for each image detected source
            #dist,ind = tree.query(np.column_stack((np.radians(detected_sources['dec']), np.radians(detected_sources['ra']))), k=1)
            tree = BallTree(np.column_stack((np.radians(adec), np.radians(ara))), leaf_size=20, metric='haversine')
            #returned is the nearest assoc neighbor for each image detected source
            dist,ind = tree.query(np.column_stack((np.radians(ddec), np.radians(dra))), k=1)
            

        cur.close()


        # Track assoc_id of matches
        matched_assoc_id=[]
        
        # Create dictionary of src_id keys & associated values
        rowdict = {}
        for i in range(len(ind)):
            sep = 3600*np.degrees(dist[i][0]) #arcsec
            match=False
            if sep < imobj.bmaj:
                match=True
                matched_assoc_id.append(assoc_sources[ind[i][0]].id)
            rowdict[detected_sources[i].src_id] = [assoc_sources[ind[i][0]].id,sep,match]
                
        # Find the unmatched associated sources
        for asrc in assoc_sources:
            if asrc.id not in matched_assoc_id:
                assoc_unmatched.append(asrc)

        # Check the detected sources
        for src in detected_sources:
            # Get the associated source object
            asrc = [msrc for msrc in assoc_sources if \
                    msrc.id == rowdict[src.src_id][0]][0]
            if rowdict[src.src_id][2]:
                # It's a match!
                src.assoc_id = asrc.id
                detected_matched.append(src)
                # Compute weighted averages
                cur_sigra_sq = asrc.e_ra * asrc.e_ra
                cur_sigdec_sq = asrc.e_dec * asrc.e_dec
                asrc.e_ra = np.sqrt(1. / (
                    (1. / cur_sigra_sq) + (1. / (src.e_ra * src.e_ra))))
                asrc.ra = (asrc.e_ra * asrc.e_ra) * (
                    (asrc.ra / cur_sigra_sq) + (src.ra / (
                        src.e_ra * src.e_ra)))
                asrc.e_dec = np.sqrt(1. / (
                    (1. / cur_sigdec_sq) + (1. / (src.e_dec * src.e_dec))))
                asrc.dec = (asrc.e_dec * asrc.e_dec) * (
                    (asrc.dec / cur_sigdec_sq) + (src.dec / (
                        src.e_dec * src.e_dec)))
                #now fluxes, they were PBCORed in Stage 2
                cur_sigtotal_sq = asrc.e_ave_total * asrc.e_ave_total
                cur_sigpeak_sq  = asrc.e_ave_peak * asrc.e_ave_peak
                asrc.e_ave_total = np.sqrt(1. / (
                    (1. / cur_sigtotal_sq) + (1. / (src.e_total_flux * src.e_total_flux))))
                asrc.e_ave_peak  = np.sqrt(1. / (
                    (1. / cur_sigpeak_sq) + (1. / (src.e_peak_flux * src.e_peak_flux))))
                asrc.ave_total = (asrc.e_ave_total * asrc.e_ave_total) * (
                    (asrc.ave_total / cur_sigtotal_sq) + (src.total_flux / (
                        src.e_total_flux * src.e_total_flux)))
                asrc.ave_peak  = (asrc.e_ave_peak * asrc.e_ave_peak) * (
                    (asrc.ave_peak / cur_sigpeak_sq) + (src.peak_flux / (
                        src.e_peak_flux * src.e_peak_flux)))
                ###
                asrc.ndetect += 1
                #one of these src.nx should be =1 the other two =0:
                asrc.ns += src.ns
                asrc.nc += src.nc
                asrc.nm += src.nm
                assoc_matched.append(asrc)
            else:
                # No match -- new source
                src.res_class = res_class
                src.ndetect = 1
                detected_unmatched.append(src)
            if not save:
                match_logger.info('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    src.src_id, rowdict[src.src_id][2], asrc.id, asrc.ra,
                    asrc.e_ra, asrc.dec, asrc.e_dec, rowdict[src.src_id][1],
                    asrc.ndetect))

    match_logger.info(' -- number of matches: {}'.format(len(detected_matched)))
    match_logger.info(' -- number of new sources to add: {}'.format(
        len(detected_unmatched)))

    return detected_matched, detected_unmatched, assoc_matched, assoc_unmatched


def filter_catalogs(conn, catalogs, imobj):
    """Selects only radio catalogs with a spatial resolution that
    lies in similar range as the current image's resolution.
    The A & B configuration equivalent resolution ranges are
    combined into one big range to include more all-sky survey catalogs.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    catalogs : list
        List of catalog names to check.
    imobj : ``database.dbclasses.Image`` instance
        Image object whose attributes are used for setting
        spatial resolution range

    Returns
    -------
    filtered_catalogs : list
        Names of catalogs which have resolutions adequate
        to proceed with the positional cross-matching.
    """
    if imobj.vcss or imobj.mosaic:
        use_range = (15.,30.) #VCSS snapshots & mosaics, regardless of config
    elif imobj.config == 'A':
        use_range = (0.,26.) #to include TGSS & LOTSS
    elif imobj.config == 'B':
        use_range = (15.,30.)
    elif imobj.config == 'C':
        use_range = (39.,80.)
    elif imobj.config == 'D':
        use_range = (99.,9999.)
    else:
        match_logger.info('Filter_catalogs received config = {} but '
                          'only A, B, C, or D are allowed...'.format(imobj.config))
        use_range = (0,0)
    ###############################

    # Find all catalogs that fall into the adequate resolution range
    cur = conn.cursor()
    filtered_catalogs = []
    for catalog in catalogs:
        try:
            catalog_res = catalogio.catalog_dict[catalog]['resolution']
        except KeyError:
            cur.execute('''SELECT resolution FROM radcat.catalogs
                WHERE name = %s''', (catalog, ))
            catalog_res = cur.fetchone()[0]
        if use_range[0] < catalog_res <= use_range[1]:
            filtered_catalogs.append(catalog)

    cur.close()

    return filtered_catalogs


def catalogmatch(conn, sources, catalog, imobj, search_radius, save):
    """Matches VLITE sources to sources from other radio
    sky survey catalogs.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        List of DetectedSource objects which need a
        sky survey catalog match.
    catalog : str
        Name of the sky survey catalog whose sources
        are being used for cross-matching.
    imobj : ``database.dbclasses.Image`` instance
        Image object whose attributes are used for setting
        the cone search center.  Only used if match_in_db
        is ``False``.
    search_radius : float
        Size of the circular search region in degrees. Only
        used if match_in_db is ``False``.
    save : bool
        If ``True``, store the catalog cross-matching results
        for future insertion into the database. If ``False``,
        print the results to the console and/or log file.
        This parameter is automatically set to ``False``
        if only the source finding and catalog matching
        stages are turned on.

    Returns
    -------
    sources : list
        DetectedSource objects with updated 'nmatches' attribute.
    catalog_matched : list
        CatalogSource objects which have been successfully
        matched to the VLITE sources.
    """
    catalog_matched = []

    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    
    match_logger.info('Attempting to match {} sources from this image to '
                      'the {} sky catalog...'.format(len(sources), catalog))

    # Print results without saving to database
    if not save:
        # Dump sources into a temporary table
        sql = (
            '''
            CREATE TEMP TABLE temp_source (
                src_id INTEGER,
                ra DOUBLE PRECISION,
                dec DOUBLE PRECISION
            );
            ''')
        cur.execute(sql)
        conn.commit()
        for src in sources:
            cur.execute('''INSERT INTO temp_source (
                src_id, ra, dec) VALUES (%s, %s, %s)''', (
                    src.src_id, src.ra, src.dec))
        conn.commit()
        # Find nearest neighbor within FOV & "match" if within half a beam
        sql = '''SELECT a.src_id, bb.id AS catalog_src_id,
            3600*q3c_dist(a.ra, a.dec, bb.ra, bb.dec) AS sep,
            3600*q3c_dist(a.ra, a.dec, bb.ra, bb.dec) < %s AS match
            FROM temp_source AS a, LATERAL (
            SELECT b.* FROM radcat.{} AS b
            WHERE q3c_join(a.ra, a.dec, b.ra, b.dec, %s)
            ORDER BY q3c_dist(a.ra, a.dec, b.ra, b.dec) ASC LIMIT 1) AS bb'''
        values = (0.5*imobj.bmin, 2.*imobj.radius)
        cur.execute(psycopg2.sql.SQL(sql).format(
            psycopg2.sql.Identifier(catalog)), values)
        rows = cur.fetchall()
        cur.execute('DROP TABLE temp_source')
        conn.commit()

        match_logger.info('-------------------------------------------------'
                          '-------------')
        match_logger.info('VLITE_src_id match catalog_src_id '
                          'separation (arcsec)')
        match_logger.info('-------------------------------------------------'
                          '-------------') 
        for row in rows:
            if row['match']:
                catalog_matched.append(row['catalog_src_id'])
            match_logger.info('{}\t\t{}\t{}\t{}'.format(
                row['src_id'], row['match'], row['catalog_src_id'], row['sep']))

    # Store results for insertion into database
    else:
        # Skip the sources which already have results for this catalog
        # (from a different image)
        assoc_ids = []
        for src in sources:
            already_matched = dbio.check_catalog_match(conn, src.id, catalog)
            if already_matched:
                continue
            else:
                assoc_ids.append(src.id)
        match_logger.info(' -- found previous matching results for {} sources'.
                          format(len(sources) - len(assoc_ids)))

        # Find nearest neighbor within half a beam
        sql = '''SELECT a.id AS assoc_id, bb.*, 
            3600*q3c_dist(a.ra, a.dec, bb.ra, bb.dec) AS sep
            FROM assoc_source AS a, LATERAL (
            SELECT b.* FROM radcat.{} AS b
            WHERE a.id IN %s AND q3c_join(a.ra, a.dec, b.ra, b.dec, %s)
            ORDER BY q3c_dist(a.ra, a.dec, b.ra, b.dec) ASC LIMIT 1) AS bb'''
        values = (tuple(assoc_ids), (0.5*(imobj.bmin/3600.)))
        cur.execute(psycopg2.sql.SQL(sql).format(
            psycopg2.sql.Identifier(catalog)), values)
        rows = cur.fetchall()

        matched_ids = []
        for row in rows:
            matched_ids.append(row['assoc_id'])
            csrc = catalogio.CatalogSource()
            dbclasses.dict2attr(csrc, row)
            catalog_matched.append(csrc)

        for src in sources:
            if src.id in matched_ids:
                # Found a match!
                try:
                    src.nmatches += 1
                except TypeError:
                    src.nmatches = 1
            else:
                if src.nmatches is None:
                    src.nmatches = 0

    cur.close()

    match_logger.info (' -- number of matches: {}'.format(len(catalog_matched)))

    return sources, catalog_matched


def check_clean(conn, sources, imobj):
    """Calculates which sources were CLEANed.
    If src within half BMIN of a CLEAN
    component "clean" is assigned True
    
    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    sources : list
        List of DetectedSource objects which need 
        clean check.
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with attribute values
        set from header info.

    Returns
    -------
    imobj : ``database.dbclasses.Image`` instance
        Initialized Image object with updated attributes
        from check clean results.
    sources : list
        List of ``database.dbclasses.DetectedSource`` objects
        updated with check clean results
    """
    #if no clean components then do nothing
    if imobj.cc is None:
        match_logger.info(' -- no CLEAN components, returning.')
        return imobj, sources
    
    cur = conn.cursor()
    
    # Dump CCs into temporary table
    sql = (
        '''
        CREATE TEMP TABLE temp_source (
        ra DOUBLE PRECISION,
        dec DOUBLE PRECISION
        );
        ''')
    cur.execute(sql)
    conn.commit()
    for cc in imobj.cc:
        cur.execute('''INSERT INTO temp_source (
            ra, dec) VALUES (%s, %s)''', (cc[0], cc[1]))
    conn.commit()
    # Find CC within half a beam of each source
    ncln=0
    match = 0.5*imobj.bmin/3600 #deg
    for src in sources:
        sql = '''SELECT COUNT(1) FROM temp_source WHERE 
            q3c_join(ra, dec, %s, %s, %s)'''
        values = (src.ra, src.dec, match)
        cur.execute(sql, values)
        row = cur.fetchall()
        count = int(row[0][0])
        if count > 0:
            src.clean = True
            ncln+=1
        else:
            src.clean = False
    imobj.nclean = ncln
    cur.execute('DROP TABLE temp_source')
    conn.commit()

    cur.close()

    match_logger.info(' -- found {}/{} sources were CLEANed'.
                          format(ncln,imobj.nsrc))
         
    return imobj, sources



def nearestneigh(sources):
    """Calculates nearest image source for all sources in image
    Parameters
    ----------
    sources : list
        List of DetectedSource objects (sources in image)

    Returns
    -------
    sources : list
        List of ``database.dbclasses.DetectedSource`` objects
        updated nn_dist & nn_src_id
    """
    #if no neighbors then return
    if len(sources) < 2:
        return sources
    
    sra = []
    sdec = []
    for src in sources:
        sra.append(src.ra)
        sdec.append(src.dec)
    sra=np.array(sra)
    sdec=np.array(sdec)
    # Make Ball Tree out of image sources with haversine distance as the metric
    tree = BallTree(np.column_stack((np.radians(sdec), np.radians(sra))), leaf_size=20, metric='haversine')

    # Return nearest 2 neighbors for each image src, 2 because 1st is the source itself
    # dist is in radians!
    dist,ind = tree.query(np.column_stack((np.radians(sdec), np.radians(sra))), k=2)

    for i in range(len(ind)):
        sources[i].nn_src_id = sources[ind[i][1]].src_id
        sources[i].nn_dist   = 3600*np.degrees(dist[i][1]) #arcsec
        #print('nearest neighbors: ',sources[i].src_id,sources[i].nn_src_id,sources[i].nn_dist)

    return sources
