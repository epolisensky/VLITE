"""radioxmatch.py contains all the machinery for cross-matching
radio point sources.

Adapted from EP's VSLOW.py.

"""
import re
import numpy as np
import psycopg2
import psycopg2.extras
from psycopg2 import sql
from skycatalog import catalogio
from database import dbclasses, dbio
import matchfuncs


res_dict = {'A' : (0., 15.), 'B' : (15., 35.),
            'C' : (35., 60.), 'D' : (60., 9999.)}


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
    res : float
        Spatial resolution of the current image in arcseconds.

    Returns
    -------
    keep : list
        List of psycopg2 row dictionary objects after
        applying the spatial resolution filtering.
    """
    keep = []
    print('Limiting to sources in resolution class {} ({}" < BMIN <= {}")'.
          format(res_class, res_dict[res_class][0], res_dict[res_class][1]))
    for row in rows:
        if row['res_class'] == res_class:
            keep.append(row)

    print(' -- {} sources remaining'.format(len(keep)))

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


def associate(conn, detected_sources, imobj, search_radius, match_in_db):
    """Associates new sources with old sources if the center
    positions of the two sources are separated by an angular
    distance less than half the size of semi-minor axis of
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
    # Find image resolution class
    for config, res_range in res_dict.items():
        if res_range[0] < imobj.bmin <= res_range[1]:
            res_class = config
    
    # Extract all previously detected sources in the same FOV
    assoc_rows = cone_search(conn, 'assoc_source', imobj.obs_ra,
                             imobj.obs_dec, search_radius)
    print('\nExtracted {} sources from assoc_source table within {} degrees.'.
          format(len(assoc_rows), search_radius))
    # Limit to sources taken from images of similar resolution
    if len(assoc_rows) > 0:
        filtered_assoc_rows = filter_res(assoc_rows, res_class)
    else:
        filtered_assoc_rows = []

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
        for asrc in filtered_assoc_rows:
            assoc_ids.append(asrc['id'])
            assoc_sources.append(dbclasses.DetectedSource())
            dbclasses.dict2attr(assoc_sources[-1], asrc)
        print('Attempting to match {} sources from this image to {} sources '
              'previously detected in VLITE images...'.format(
                  len(detected_sources), len(assoc_sources)))

        detected_matched = []
        detected_unmatched = []
        assoc_matched = []
        assoc_unmatched = []

        if not match_in_db:
            print('-----------------------------------------------------------'
                  '--------------------------------------------------------')
            print('src_id match assoc_id\tra\t\te_ra\t\t\tdec\t\te_dec\t\t'
                  'separation (arcsec)')
            print('-----------------------------------------------------------'
                  '--------------------------------------------------------')
            for src in detected_sources:
                match, asrc, min_sep = matchfuncs.simple_match(
                    src, assoc_sources, imobj.bmin)
                if match:
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
                    asrc.ndetect += 1
                    assoc_matched.append(asrc)
                else:
                    src.ndetect = 1
                    src.res_class = res_class
                    detected_unmatched.append(src)
                try:
                    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        src.src_id, match, asrc.id, asrc.ra, asrc.e_ra,
                        asrc.dec, asrc.e_dec, min_sep))
                except AttributeError:
                    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        src.src_id, match, None, src.ra, src.e_ra,
                        src.dec, src.e_dec, min_sep))
        else:
            # Match sources if they are separatred by less than half a beam
            cur = conn.cursor()
            #sql = '''SELECT a.src_id, b.id FROM detected_source AS a, 
            #    assoc_source AS b WHERE a.image_id = %s AND b.id IN %s AND
            #    q3c_join(a.ra, a.dec, b.ra, b.dec, %s)'''
            sql = '''SELECT a.src_id, bb.id FROM detected_source AS a,
                LATERAL ( SELECT b.* FROM assoc_source AS b
                WHERE a.image_id = %s AND b.id IN %s AND
                q3c_join(a.ra, a.dec, b.ra, b.dec, %s) ORDER BY
                q3c_dist(a.ra, a.dec, b.ra, b.dec) ASC LIMIT 1) AS bb'''
            values = (imobj.id, tuple(assoc_ids), (0.5*(imobj.bmin/3600.)))
            cur.execute(sql, values)
            rows = cur.fetchall()
            cur.close()

            # Create dictionary of src_id keys and assoc_id values of matches
            rowdict = {}
            for row in rows:
                rowdict[row[0]] = row[1]

            assoc_unmatched = [asrc for asrc in assoc_sources if \
                               asrc.id not in rowdict.values()]

            for src in detected_sources:
                if src.src_id in rowdict.keys():
                    # It's a match!
                    src.assoc_id = rowdict[src.src_id]
                    detected_matched.append(src)
                    # Update the associated source
                    asrc = [msrc for msrc in assoc_sources if \
                            msrc.id == src.assoc_id][0]
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
                    asrc.ndetect += 1
                    assoc_matched.append(asrc)
                else:
                    # No match -- new source
                    src.res_class = res_class
                    src.ndetect = 1
                    detected_unmatched.append(src)

    print(' -- number of matches: {}'.format(len(detected_matched)))
    print(' -- number of new sources to add: {}'.format(
        len(detected_unmatched)))

    return detected_matched, detected_unmatched, assoc_matched, assoc_unmatched


def filter_catalogs(conn, catalogs, res):
    """Selects only radio catalogs with a spatial resolution that
    lies in the same range as the current image's resolution.
    The A & B configuration equivalent resolution ranges are
    combined into one big range to include more all-sky survey catalogs.

    Parameters
    ----------
    conn : ``psycopg2.extensions.connect`` instance
        The PostgreSQL database connection object.
    catalogs : list
        List of catalog names to check.
    res : float
        Spatial resolution of the current image.

    Returns
    -------
    filtered_catalogs : list
        Names of catalogs which have resolutions adequate
        to proceed with the positional cross-matching.
    """   
    # Determine which resolution range the image belongs in
    for config, res_range in res_dict.items():
        if res_range[0] < res <= res_range[1]:
            use_range = res_range
            # Combine highest resolutions to allow for more catalogs
            if config == 'A' or config == 'B':
                use_range = (res_dict['A'][0], res_dict['B'][1])

    # Find all catalogs that fall into the adequate resolution range
    cur = conn.cursor()
    filtered_catalogs = []
    for catalog in catalogs:
        try:
            catalog_res = catalogio.catalog_dict[catalog]['resolution']
        except KeyError:
            cur.execute('''SELECT resolution FROM skycat.catalogs
                WHERE name = %s''', (catalog, ))
            catalog_res = cur.fetchone()[0]
        if use_range[0] < catalog_res <= use_range[1]:
            filtered_catalogs.append(catalog)

    cur.close()

    return filtered_catalogs


def catalogmatch(conn, sources, catalog, imobj, search_radius, match_in_db):
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
    match_in_db : bool
        If ``True``, do the cross-matching inside the database
        between the **assoc_source** and catalog tables. If ``False``,
        do the cross-matching outside the database by first extracting
        the catalog sources using a cone search. This option should
        be set to ``False`` when the VLITE sources to be matched are
        not in the **assoc_source** table, as is the case when results
        are not saved to the database.

    Returns
    -------
    sources : list
        DetectedSource objects with updated 'nmatches' attribute.
    catalog_matched : list
        CatalogSource objects which have been successfully
        matched to the VLITE sources.
    """
    catalog_matched = []
    
    # Cross-match outside the DB
    if not match_in_db:
        # Extract catalog sources
        catalog_rows = cone_search(conn, catalog, imobj.obs_ra, imobj.obs_dec,
                                   search_radius, 'skycat')
        print('\nExtracted {} sources from {} within {} degrees.'.
              format(len(catalog_rows), catalog, search_radius))
        # Sky survey does not cover this sky region, move on
        if not catalog_rows:
            return
        
        # Translate row dictionaries to CatalogSource objects
        catalog_sources = []
        for csrc in catalog_rows:
            catalog_sources.append(catalogio.CatalogSource())
            dbclasses.dict2attr(catalog_sources[-1], csrc)
        print('Attempting to match {} sources from this image to {} '
              'sources from the {} sky catalog...'.format(
                  len(sources), len(catalog_sources), catalog))
        
        print('-----------------------------------------------------')
        print('VLITE_src_id match catalog_src_id separation (arcsec)')
        print('-----------------------------------------------------')
        
        for src in sources:
            match, csrc, min_sep = matchfuncs.simple_match(
                src, catalog_sources, imobj.bmin)
            if match: # Found a match!
                try:
                    src.nmatches += 1
                except TypeError:
                    src.nmatches = 1
                csrc.assoc_id = src.id
                csrc.sep = min_sep
                catalog_matched.append(csrc)
            else:
                if src.nmatches is None:
                    src.nmatches = 0
            print('{}\t{}\t{}\t{}'.format(src.src_id, match, csrc.id, csrc.sep))


    # Cross-match inside the DB
    else:
        print('\nAttempting to match {} sources from this image to '
              'the {} sky catalog...'.format(len(sources), catalog))
        # Skip the sources which already have results for this catalog
        # (from a different image)
        assoc_ids = []
        for src in sources:
            already_matched = dbio.check_catalog_match(conn, src.id, catalog)
            if already_matched:
                continue
            else:
                assoc_ids.append(src.id)
        print(' -- found previous matching results for {} sources'.format(
            len(sources) - len(assoc_ids)))

        # Find nearest neighbor within half a beam
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        sql = '''SELECT a.id AS assoc_id, bb.*,
            3600*q3c_dist(a.ra, a.dec, bb.ra, bb.dec) AS sep
            FROM assoc_source AS a, LATERAL (
            SELECT b.* FROM skycat.{} AS b
            WHERE a.id IN %s AND q3c_join(a.ra, a.dec, b.ra, b.dec, %s)
            ORDER BY q3c_dist(a.ra, a.dec, b.ra, b.dec) ASC LIMIT 1)
            AS bb'''
        values = (tuple(assoc_ids), (0.5*(imobj.bmin/3600.)))
        cur.execute(psycopg2.sql.SQL(sql).format(
            psycopg2.sql.Identifier(catalog)), values)
        rows = cur.fetchall()
        cur.close()

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

    print(' -- number of matches: {}'.format(len(catalog_matched)))

    return sources, catalog_matched
