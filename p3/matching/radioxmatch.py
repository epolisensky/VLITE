"""radioXmatch.py contains all the machinery for cross-matching
sources detected in images to other radio frequency catalogs.
It currently supports source catalog input as a list of
DetectedSource objects (from pybdsf_source.py) or as database.
Radio survey catalogs are read in from the SkyCatalogs
database. Lists of the sources with matches, sources without
matches, and catalog sources what were matched are returned.

Adapted from EP's VSLOW.py.

Post-Processing Pipeline (P3) Stage 4"""


import numpy as np
import psycopg2
import psycopg2.extras
from psycopg2 import sql
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from skycatalog import catalogio
from sourcefinding import beam_correction
from database import dbclasses
import matchfuncs


def write_region(srclist, impath, ext='.reg'):
    """Writes a ds9 regions file from given source list."""
    fname = impath[:-5] + ext
    with open(fname, 'w') as f:
        f.write('global color=red font="helvetica 10 normal" '
                'select=1 highlite=1 edit=1 move=1 delete=1 '
                'include=1 fixed=0 source\n')
        f.write('fk5\n')
        for src in srclist:
            f.write('ellipse(%f,%f,%.2f",%.2f",%.1f) # text={%s}\n' % (
                src.ra, src.dec, src.maj, src.min, src.pa + 90.0, src.name))


def search_radius(imobj):
    imsize = float(imobj.imsize.split()[1].split(')')[0])
    fov = (imsize * imobj.pixel_scale) / 3600. # deg
    return fov / 2.


def estimate_sn(conn, assoc_sources, imobj):
    imdata, hdr = imobj.read()
    data = imdata.reshape(imdata.shape[2:])
    wcs = WCS(hdr).celestial
    dxdy = int(round(imobj.rms_box[0] / 2.))

    cur = conn.cursor()

    for asrc in assoc_sources:
        center = SkyCoord(imobj.obs_ra, imobj.obs_dec, unit='deg')
        srcloc = SkyCoord(asrc.ra, asrc.dec, unit='deg')
        angdist = center.separation(srcloc)
        # Correct for beam response - only 1D (sym. beam) for now
        pbcorr = beam_correction.find_nearest_pbcorr(angdist.degree)
        # Estimate noise
        pixcoords = wcs.wcs_world2pix([[asrc.ra, asrc.dec]], 1)
        x = int(round(pixcoords[0,0])) # column
        y = int(round(pixcoords[0,1])) # row
        noise = np.std(data[y-dxdy:y+dxdy, x-dxdy:x+dxdy]) / pbcorr
        # Calculate median peak flux
        cur.execute('SELECT peak_flux FROM raw_source WHERE assoc_id = %s',
                    (asrc.id, ))
        pkfluxes = cur.fetchall()
        # PB correct for now -- in the future this will be a separate col
        peak_flux = np.median(pkfluxes) / pbcorr
        asrc.sn = peak_flux / noise

    cur.close()
    return assoc_sources


def limit_res(rows, res, res_tol):
    print('\nLimiting to sources with resolution = '
          '{:.1f}+/-{}"'.format(res, res_tol))
    keep = []
    for row in rows:
        if res - res_tol < row['beam'] < res + res_tol:
            keep.append(row)

    print(' -- {} sources remaining'.format(len(keep)))

    return keep


def check_previous(conn, src, search_radius, res_tol):
    resup, reslo = src.beam + res_tol, src.beam - res_tol
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute('''SELECT id FROM image 
        WHERE q3c_join(%s, %s, obs_ra, obs_dec, %s) AND
        bmaj BETWEEN %s AND %s''',
                (src.ra, src.dec, search_radius, reslo, resup))
    prev_images = cur.fetchall()
    cur.close()

    return prev_images


def cone_search(conn, schema, table, center_ra, center_dec, radius):
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
    

def box_query(conn, schema, table, box):
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    query = '''SELECT * FROM {}.{} WHERE q3c_poly_query(ra, dec, 
        '{%s, %s, %s, %s, %s, %s, %s, %s}')'''
    cur.execute(psycopg2.sql.SQL(query).format(
        psycopg2.sql.Identifier(schema), psycopg2.sql.Identifier(table)),
                (box[0,0], box[0,1],
                 box[1,0], box[1,1],
                 box[2,0], box[2,1],
                 box[3,0], box[3,1]))
    rows = cur.fetchall()
    cur.close()

    return rows


def associate(conn, detected_sources, imobj, search_radius, res_tol):
    """Associates sources.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    box : 4x2 numpy array
        Array containing RA, Dec of each corner of the box defining
        the area of the image used for source finding. Corners are
        listed starting at the bottom left corner and go clockwise.
    detected_source : list of `DetectedSource` objects
        List of the sources extracted from the current image.
    imobj : database.dbclasses.Image instance
        Initialized `Image` object with attribute values
        set from header info & updated with source finding results.
    res_tol : float
        Defines the acceptable range of spatial resolutions in arcsec
        when considering a source from the assoc_source table for
        cross-matching.

    Returns
    -------
    detected_matched : list of `DetectedSource` objects
        Sources extracted from the image that were successfully
        associated with previously detected sources stored in the
        assoc_source table. The assoc_id column for each of these
        sources is updated in the detected_source table with the id of
        the corresponding row/source in the assoc_source table.
    detected_unmatched : list of `DetectedSource` objects
        Sources extracted from the image which could not be successfully
        associated with previously detected sources. These are added
        to the assoc_table as new detections and then cross-matched
        with other sky survey catalogs. If no match is found, then
        the source becomes a transient candidate.
    assoc_matched : list of `DetectedSource` objects
        Previously detected sources from the assoc_source table which
        were successfully associated with sources from the new image. 
        Size and position properties of these sources are updated in the
        assoc_source table to reflect the weighted average of all 
        detections. The number of detections (ndetect) is increased 
        by one.
    assoc_unmatched : list of `DetectedSource` objects
        Previously detected sources from the assoc_source table which
        were not successfully associated with sources from the new image.
        Right now, this list is being returned for posterity and is not
        actually used for anything.
    imobj : database.dbclasses.Image instance
        Initialized `Image` object with updated stage attribute.
    """
    #assoc_rows = box_query(conn, 'public', 'assoc_source', box)
    assoc_rows = cone_search(conn, 'public', 'assoc_source', imobj.obs_ra,
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
            if match == 1: # Found a match!
                src.assoc_id = asrc.id
                detected_matched.append(src)
                # Compute (weighted, eventually) averages
                attrs = ['ra', 'e_ra', 'dec', 'e_dec', 'maj', 'e_maj',
                         'min', 'e_min', 'pa', 'e_pa']
                for attr in attrs:
                    prevattr = getattr(asrc, attr)
                    curattr = getattr(src, attr)
                    setattr(asrc, attr, ((prevattr * asrc.ndetect) \
                                         + curattr) / (asrc.ndetect + 1))
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
    # Extract catalog sources
    catalog_rows = cone_search(conn, 'skycat', catalog, imobj.obs_ra,
                               imobj.obs_dec, search_radius)
    #table = 'skycat.' + catalog
    #catalog_rows = box_query(conn, table, box)
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
            if match == 1: # Found a match!
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
