"""radioxmatch.py contains all the machinery for cross-matching
radio point sources.

Adapted from EP's VSLOW.py.

Post-Processing Pipeline (P3) Stage 4

"""
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
        (i.e. '1.5GHz.0137+331.IPln1.matches.reg'). Default
        is '.reg'.
    """
    fname = impath[:-5] + ext
    with open(fname, 'w') as f:
        f.write('global color=red font="helvetica 10 normal" '
                'select=1 highlite=1 edit=1 move=1 delete=1 '
                'include=1 fixed=0 source\n')
        f.write('fk5\n')
        for src in srclist:
            f.write('ellipse(%f,%f,%.2f",%.2f",%.1f) # text={%s}\n' % (
                src.ra, src.dec, src.maj, src.min, src.pa + 90.0, src.name))


def limit_res(rows, res, res_tol):
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
    res_tol : float
        Acceptable range above and below the given resolution
        for sources to be considered for matching.

    Returns
    -------
    keep : list
        List of `psycopg2` row dictionary objects after
        applying the spatial resolution filtering.
    """
    print('\nLimiting to sources with resolution = '
          '{:.1f}+/-{}"'.format(res, res_tol))
    keep = []
    for row in rows:
        if res - res_tol < row['beam'] < res + res_tol:
            keep.append(row)

    print(' -- {} sources remaining'.format(len(keep)))

    return keep


def check_previous(conn, src, search_radius, res_tol):
    """Searches the database image table for images
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
    res_tol : float
        Spatial resolution tolerance to define range
        of acceptable resolutions to consider.

    Returns
    -------
    prev_images : list
        List of ids of the images which could have
        contained the given point.
    """
    resup, reslo = src.beam + res_tol, src.beam - res_tol
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute('''SELECT id FROM image 
        WHERE q3c_join(%s, %s, obs_ra, obs_dec, %s) AND
        bmaj BETWEEN %s AND %s''',
                (src.ra, src.dec, search_radius, reslo, resup))
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
    schema : str, optional
        The database schema which contains the table.
        Default is 'public'.
    table : str
        Name of the database table to query.
    center_ra : float
        Right ascension coordinate of the circle center
        in degrees.
    center_dec : float
        Declination coordinate of the circle center in degrees.
    radius : float
        Size of the circular search region in degrees.

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
    

def associate(conn, detected_sources, imobj, search_radius, res_tol):
    """Associates new sources with old sources.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    detected_sources : list of DetectedSource objects
        List of the sources extracted from the current image.
    imobj : database.dbclasses.Image instance
        Initialized Image object with attribute values
        set from header info & updated with source finding results.
    search_radius : float
        Size of the circular search region in degrees.
    res_tol : float
        Defines the acceptable range of spatial resolutions in arcsec
        when considering a source from the assoc_source table for
        cross-matching.

    Returns
    -------
    detected_matched : list of DetectedSource objects
        Sources extracted from the image that were successfully
        associated with previously detected sources stored in the
        assoc_source table.
    detected_unmatched : list of DetectedSource objects
        Sources extracted from the image which could not be successfully
        associated with previously detected sources. These are added
        to the assoc_table as new detections and then cross-matched
        with other sky survey catalogs.
    assoc_matched : list of DetectedSource objects
        Previously detected sources from the assoc_source table which
        were successfully associated with sources from the new image. 
        Size and position properties of these sources are updated in the
        assoc_source table to reflect the weighted average of all 
        detections. The number of detections (ndetect) is increased 
        by one.
    assoc_unmatched : list of DetectedSource objects
        Previously detected sources from the assoc_source table which
        were not successfully associated with sources from the new image.
        Any of these non-detections with no sky catalog matches
        (nmatches = 0) are recorded in the vlite_unique table.
    """
    assoc_rows = cone_search(conn, 'assoc_source', imobj.obs_ra,
                             imobj.obs_dec, search_radius)
    print('\nExtracted {} sources from assoc_source table within {} degrees.'.
          format(len(assoc_rows), search_radius))
    if not assoc_rows:
        # No previous sources found in that sky region
        for src in detected_sources:
            src.beam = imobj.bmaj
            src.ndetect = 1
        detected_matched = []
        detected_unmatched = detected_sources
        assoc_matched = []
        assoc_unmatched = []
    else:
        # Limit to sources taken from images of similar resolution
        limited_assoc_rows = limit_res(assoc_rows, imobj.bmaj, res_tol)
        # Translate row dictionaries to DetectedSource objects
        assoc_sources = []
        for asrc in limited_assoc_rows:
            assoc_sources.append(dbclasses.DetectedSource())
            dbclasses.dict2attr(assoc_sources[-1], asrc)
        print('\nAttempting to match {} sources from this image to {} sources '
              'previously detected in VLITE images...'.format(
                  len(detected_sources), len(assoc_sources)))
        detected_matched = []
        detected_unmatched = []
        assoc_matched = []
        assoc_unmatched = []
        for src in detected_sources:
            match, asrc, min_der = matchfuncs.deRuitermatch(
                src, assoc_sources, imobj.bmaj)
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
                src.beam = imobj.bmaj
                src.ndetect = 1
                detected_unmatched.append(src)

        # asrc is not returned when there is no match, so loop outside
        for asrc in assoc_sources:
            if asrc not in assoc_matched:
                assoc_unmatched.append(asrc)

    print(' -- number of matches: {}'.format(len(detected_matched)))
    print(' -- number of new sources to add: {}'.format(
        len(detected_unmatched)))
    print(' -- number of unmatched previously detected sources: {}'.format(
        len(assoc_unmatched)))

    return detected_matched, detected_unmatched, assoc_matched, assoc_unmatched
            

def catalogmatch(conn, sources, catalog, imobj, search_radius):
    """Matches VLITE sources to other source from
    other sky surveys.

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
        print('\nAttempting to match {} sources from this image to {} '
              'sources from the {} sky catalog...'.format(
                  len(sources), len(catalog_sources), catalog))
        catalog_matched = []
        sources_unmatched = []
        for src in sources:
            match, csrc, min_der = matchfuncs.deRuitermatch(
                src, catalog_sources, imobj.bmaj)
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
