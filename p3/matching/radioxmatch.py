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


def limit_res(rows, res):
    """Filters out sources extracted from the database
    which originate from images with spatial resolutions
    outside the acceptable range.

    Parameters
    ----------
    rows : list
        `psycopg2` row dictionary objects extracted
        from the database.
    res : float
        Spatial resolution of the current image in arcseconds.

    Returns
    -------
    keep : list
        List of `psycopg2` row dictionary objects after
        applying the spatial resolution filtering.
    """
    res = int(round(res))
    keep = []
    # A/B+
    if res <= 15.0:
        print('Limiting to sources with BMIN <= 15"')
        for row in rows:
            if row['beam'] <= 15.0:
                keep.append(row)
    # B
    elif 15. < res <= 35.:
        print('Limiting to sources with 15" < BMIN <= 35"')
        for row in rows:
            if 15. < row['beam'] <= 35.:
                keep.append(row)
    # C
    elif 35. < res <= 60.:
        print('Limiting to sources with 35" < BMIN <= 60"')
        for row in rows:
            if 35. < row['beam'] <= 60.:
                keep.append(row)
    # D
    else:
        print('Limiting to sources with BMIN > 60"')
        for row in rows:
            if row['beam'] > 60.:
                keep.append(row)

    print(' -- {} sources remaining'.format(len(keep)))

    return keep


def check_previous(conn, src, search_radius):
    """Searches the database **image** table for images
    of similar spatial resolution which cover an area
    on the sky containing a given point.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    src : database.dbclasses.DetectedSource instance
        Source whose position will be used to find
        all previously processed images which could
        have contained it.
    search_radius : float
        Radius defining the size of the circular search
        area in degrees.

    Returns
    -------
    prev_images : list
        List of ids of the images which could have
        contained the given point.
    """
    if src.beam <= 15.0:
        reslo, reshi = 0., 15.
    elif 15. < src.beam <= 35.:
        reslo, reshi = 15., 35.
    elif 35. < src.beam <= 60.:
        reslo, reshi = 35., 60.
    else:
        reslo, reshi = 60., 99999.

    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute('''SELECT id FROM image 
        WHERE q3c_join(%s, %s, obs_ra, obs_dec, %s) AND
        bmin > %s AND bmin <= %s''',
                (src.ra, src.dec, search_radius, reslo, reshi))
    prev_images = cur.fetchall()
    cur.close()

    return prev_images


def cone_search(conn, table, center_ra, center_dec, radius, schema='public'):
    """Extracts all sources from the specified database
    table which fall within a circular region on the sky.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
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
    

def associate(conn, detected_sources, imobj, search_radius):
    """Associates new sources with old sources.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    detected_sources : list
        List of the sources extracted from the current image.
    imobj : database.dbclasses.Image instance
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
        with other sky survey catalogs.
    assoc_matched : list
        Previously detected sources from the **assoc_source** table which
        were successfully associated with sources from the new image. 
        Positions of these sources are updated in the **assoc_source**
        table to reflect the weighted average of all detections. The 
        number of detections (ndetect) is increased by one.
    assoc_unmatched : list
        Previously detected sources from the **assoc_source** table which
        were not successfully associated with sources from the new image.
        Any of these non-detections with no sky catalog matches
        (nmatches = 0) are recorded in the **vlite_unique** table.
    """
    # Extract all previously detected sources in the same FOV
    assoc_rows = cone_search(conn, 'assoc_source', imobj.obs_ra,
                             imobj.obs_dec, search_radius)
    print('\nExtracted {} sources from assoc_source table within {} degrees.'.
          format(len(assoc_rows), search_radius))
    
    if not assoc_rows:
        # No previous sources found in that sky region
        for src in detected_sources:
            src.beam = imobj.bmin
            src.ndetect = 1
        detected_matched = []
        detected_unmatched = detected_sources
        assoc_matched = []
        assoc_unmatched = []
    else:
        # Limit to sources taken from images of similar resolution
        limited_assoc_rows = limit_res(assoc_rows, imobj.bmin)
        # Translate row dictionaries to DetectedSource objects
        assoc_sources = []
        for asrc in limited_assoc_rows:
            assoc_sources.append(dbclasses.DetectedSource())
            dbclasses.dict2attr(assoc_sources[-1], asrc)
        print('Attempting to match {} sources from this image to {} sources '
              'previously detected in VLITE images...'.format(
                  len(detected_sources), len(assoc_sources)))
        detected_matched = []
        detected_unmatched = []
        assoc_matched = []
        assoc_unmatched = []
        # Use de Ruiter radius to cross-match each new source with old ones
        for src in detected_sources:
            match, asrc, min_der = matchfuncs.deRuitermatch(
                src, assoc_sources, imobj.bmin)
            if match: # Found a match!
                src.assoc_id = asrc.id
                detected_matched.append(src)
                # Compute weighted averages
                cur_sigra_sq = asrc.e_ra * asrc.e_ra
                cur_sigdec_sq = asrc.e_dec * asrc.e_dec
                asrc.e_ra = np.sqrt(1. / (
                    (1. / cur_sigra_sq) + (1. / (src.e_ra * src.e_ra))))
                asrc.ra = (asrc.e_ra * asrc.e_ra) * (
                    (asrc.ra / cur_sigra_sq) + (src.ra / (src.e_ra * src.e_ra)))
                asrc.e_dec = np.sqrt(1. / (
                    (1. / cur_sigdec_sq) + (1. / (src.e_dec * src.e_dec))))
                asrc.dec = (asrc.e_dec * asrc.e_dec) * (
                    (asrc.dec / cur_sigdec_sq) + (
                        src.dec / (src.e_dec * src.e_dec)))
                asrc.ndetect += 1
                assoc_matched.append(asrc)
            else:
                src.beam = imobj.bmin
                src.ndetect = 1
                detected_unmatched.append(src)

        # asrc is not returned when there is no match, so loop outside
        for asrc in assoc_sources:
            if asrc not in assoc_matched:
                assoc_unmatched.append(asrc)

    print(' -- number of matches: {}'.format(len(detected_matched)))
    print(' -- number of new sources to add: {}'.format(
        len(detected_unmatched)))

    return detected_matched, detected_unmatched, assoc_matched, assoc_unmatched
            

def catalogmatch(conn, sources, catalog, imobj, search_radius):
    """Matches VLITE sources to sources from other radio
    sky survey catalogs.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    sources : list
        List of DetectedSource objects which need a
        sky survey catalog match.
    catalog : str
        Name of the sky survey catalog whose sources
        are being used for cross-matching.
    imobj : database.dbclasses.Image instance
        Image object whose attributes are used for setting
        search center.
    search_radius : float
        Size of the circular search region in degrees.

    Returns
    -------
    sources : list
        DetectedSource objects with updated nmatches attribute.
    catalog_matched : list
        CatalogSource objects which have been successfully
        matched to the VLITE sources.
    """            
    # Extract catalog sources
    catalog_rows = cone_search(conn, catalog, imobj.obs_ra, imobj.obs_dec,
                               search_radius, 'skycat')
    print('\nExtracted {} sources from {} within {} degrees.'.
          format(len(catalog_rows), catalog, search_radius))
    if not catalog_rows:
        # Sky survey does not cover this sky region, move on
        return
    else:
        # Translate row dictionaries to CatalogSource objects
        catalog_sources = []
        for csrc in catalog_rows:
            catalog_sources.append(catalogio.CatalogSource())
            dbclasses.dict2attr(catalog_sources[-1], csrc)
        print('Attempting to match {} sources from this image to {} '
              'sources from the {} sky catalog...'.format(
                  len(sources), len(catalog_sources), catalog))
        # Loop through sources to cross-match
        catalog_matched = []
        sources_unmatched = []
        for src in sources:
            # Skip the sources which already have results for this catalog
            # (from a different image)
            already_matched = dbio.check_catalog_match(conn, src.id, catalog)
            if already_matched:
                continue
            match, csrc, min_der = matchfuncs.deRuitermatch(
                src, catalog_sources, imobj.bmin)
            if match: # Found a match!
                try:
                    src.nmatches += 1
                except TypeError:
                    src.nmatches = 1
                csrc.assoc_id = src.id
                csrc.min_deruiter = min_der
                catalog_matched.append(csrc)
            else:
                if src.nmatches is None:
                    src.nmatches = 0
                sources_unmatched.append(src)

    print(' -- number of matches: {}'.format(len(catalog_matched)))
    print(' -- number of unmatched VLITE sources: {}'.format(
        len(sources_unmatched)))

    return sources, catalog_matched
