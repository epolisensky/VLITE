"""radioXmatch.py contains all the machinery for cross-matching
sources detected in images to other radio frequency catalogs.
It currently supports source catalog input as a list of
DetectedSource objects (from pybdsf_source.py) or as database.
Radio survey catalogs are read in from the SkyCatalogs
database. Lists of the sources with matches, sources without
matches, and catalog sources what were matched are returned.

Adapted from EP's VSLOW.py.

Post-Processing Pipeline (P3) Stage 4"""


import os
import sys
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
                src['ra'], src['dec'], src['maj'], src['min'],
                src['pa'] + 90.0, src['name']))


def update_database(database, imname, matches, non_matches):
    """Updates the catalog_id, match_id, & min_deRuiter
    columns in the database Source table with the results
    from the catalog cross-matching."""
    print('\nAdding catalog cross-match results to database Source table.')
    cur, conn = dbconnect(database)
    cur.execute('SELECT id FROM Image WHERE filename = ?', (imname, ))
    img_id = cur.fetchone()[0]

    for src in matches + non_matches:
        cur.execute('''UPDATE AssocSource SET catalog_id = ?, match_id = ?,
            min_deRuiter = ? WHERE src_id = ? AND image_id = ?''',
                    (src['catalog_id'], src['match_id'], src['min_deRuiter'],
                     src['src_id'], img_id))
    conn.commit()
    cur.close()


def cone_search(dbname, tables, center_coords, search_radius):
    cur, conn = dbconnect(dbname, load_ext=True)
    cra = center_coords[0]
    cdec = center_coords[1]
    r = np.radians(search_radius)
    d2r = np.pi / 180.
    print('\nSelecting sources from {} within {} degrees of {:.3f}, '
          '{:.3f}'.format(tables, int(round(search_radius)), cra, cdec))
    tabdict = {}
    """
    if same_res:
        try:
            beam = kwargs['beam']
            res_tol = kwargs['res_tol']
            bmup = beam + res_tol
            bmlo = beam - res_tol
            print('\nLimiting to sources with resolution = {:.1f}+/-{}"'.format(
                beam, res_tol))
            for table in tables:
                query = '''SELECT * FROM %s WHERE beam BETWEEN ? AND ? AND
                    2 * ASIN(SQRT(SIN(((?-dec)/2)*?) * SIN(((?-dec)/2)*?) 
                    + COS(?*?) * COS(dec*?) * SIN(((?-ra)/2)*?) 
                    * SIN(((?-ra)/2)*?))) <= ?''' % table
                cur.execute(query, (bmlo, bmup, cdec, d2r, cdec, d2r, cdec,
                                    d2r, d2r, cra, d2r, cra, d2r, r))
                tabdict[table] = cur.fetchall()
        except KeyError:
            print('\nNo beam or resolution tolerance provided.')
            print('Overriding resolution constraint.')
            same_res = False
    """

    for table in tables:
        query = '''SELECT * FROM %s WHERE 
            2 * ASIN(SQRT(SIN(((?-dec)/2)*?) * SIN(((?-dec)/2)*?) 
            + COS(?*?) * COS(dec*?) * SIN(((?-ra)/2)*?) 
            * SIN(((?-ra)/2)*?))) <= ?''' % table
        cur.execute(query, (cdec, d2r, cdec, d2r, cdec, d2r, d2r, cra, d2r,
                            cra, d2r, r))
        tabdict[table] = cur.fetchall()
        print(' -- retrieved {} sources from {}\n'.format(len(tabdict[table]),
                                                         table))

    cur.close()

    return tabdict


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


def limit_res(assoc_rows, res, res_tol):
    print('\nLimiting to associated sources with resolution = '
          '{:.1f}+/-{}"'.format(res, res_tol))
    keep = []
    for asrc in assoc_rows:
        if res - res_tol < asrc['beam'] < res + res_tol:
            keep.append(asrc)

    print(' -- {} sources remaining'.format(len(keep)))

    return keep


def associate(conn, box, detected_sources, imobj, res_tol):
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
    assoc_matched : list of `AssociatedSource` objects
        Previously detected sources from the assoc_source table which
        were successfully associated with sources from the new image. 
        Size and position properties of these sources are updated in the
        assoc_source table to reflect the weighted average of all 
        detections. The number of detections (ndetect) is increased 
        by one.
    assoc_unmatched : list of `AssociatedSource` objects
        Previously detected sources from the assoc_source table which
        were not successfully associated with sources from the new image.
        Right now, this list is being returned for posterity and is not
        actually used for anything.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    query = '''SELECT * FROM assoc_source
        WHERE q3c_poly_query(ra, dec, '{%s, %s, %s, %s, %s, %s, %s, %s}')'''
    cur.execute(query, (box[0,0], box[0,1],
                        box[1,0], box[1,1],
                        box[2,0], box[2,1],
                        box[3,0], box[3,1]))
    assoc_rows = cur.fetchall()
    cur.close()
    print('\nExtracted {} sources from assoc_source table within {}.'.
          format(len(assoc_rows), box))
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
        # Translate row dictionaries to AssociatedSource objects
        assoc_sources = []
        for asrc in limited_assoc_rows:
            assoc_sources.append(dbclasses.AssociatedSource())
            dbclasses.dict2attr(assoc_sources[-1], asrc)
        print('\nAttempting to match {} sources from this image to {} sources '
              'previously detected in VLITE images...'.format(
                  len(detected_sources), len(assoc_sources)))
        detected_matched = []
        detected_unmatched = []
        assoc_matched = []
        assoc_unmatched = []
        for src in detected_sources:
            match, asrc = matchfuncs.deRuitermatch(src, assoc_sources,
                                                   imobj.bmaj)
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
        #if assoc_unmatched:
            # Updates nopp if unmatched asrc is above detection threshold
            #assoc_unmatched = estimate_sn(conn, assoc_unmatched, imobj)
        
    print(' -- number of matches: {}'.format(len(detected_matched)))
    print(' -- number of new sources to add: {}'.format(
        len(detected_unmatched)))
    print(' -- number of unmatched previously detected sources: {}'.format(
        len(assoc_unmatched)))

    return detected_matched, detected_unmatched, assoc_matched, assoc_unmatched  
            

def catalogmatch(image, sources, skycat, catalogs):
    # Connect the sky catalog database
    catcur, catconn = dbconnect(skycat)

    # Extract catalog sources
    catdict = catalog_extract(catcur, catalogs, search_range(image))

    # Define beam major axis FWHM:
    bmaj = image['bmaj'] / 3600. # deg

    # Cross-match each source to the catalog sources
    print('\nSearching for catalog cross-matches...')
    matches = []
    non_matches = []
    cat_matches = []
    for src in sources:
        idx = 0
        match, catsrc = matchfuncs.deRuitermatch(
            src, catdict[catalogs[idx]], bmaj)
        while match != 1:
            idx += 1
            if idx < len(catalogs):
                match, catsrc = matchfuncs.deRuitermatch(
                    src, catdict[catalogs[idx]], bmaj)
            else:
                # Ran out of catalogs to search -- no match found
                src['catalog_id'] = -1
                src['match_id'] = -1
                non_matches.append(src)
                break
        else:
            # Found a match!
            src['catalog_id'] = catsrc['catalog_id']
            src['match_id'] = catsrc['id']
            matches.append(src)
            cat_matches.append(catsrc)
            
    print('\nMatched {}/{} sources.'.format(len(matches), len(sources)))
    catcur.close()

    return matches, non_matches, cat_matches


def prep(catalogs=['NVSS'], **kwargs):
    """Prepares all possible inputs for catalog cross-matching.
    Input options are either a database of images and their sources
    (like the one written by database.py) or an ImageTable object
    and list of DetectedSource objects which are used to populate
    the database tables. If the input is the former, the database
    keyword must be set to True and the full path to the database 
    along with a list of image names for which sources are 
    to be cross-matched must be provided (i.e. database=True,
    file='/path/to/database.sqlite', images=['']). If the input is
    the latter, then the database keyword must be set to False and
    the ImageTable and DetectedSource objects must be provided
    (i.e. database=False, objects=(imgtbl, sources)). In both
    cases, the list of catalogs to check can be specified through
    the catalogs keyword with the default being only NVSS."""  
    # main(catalogs=catalogs, database=False, objects=(imgtbl, sources))
    # main(catalogs=catalogs, database=True, images=imglist,
    #      file='/data3/erichards/codes/p3/test/test.sqlite')
    if kwargs['database'] is True:
        try:
            dbname = kwargs['file']
            if not os.path.exists(dbname):
                print('ERROR: File or path {} does not exist!'.format(dbname))
                sys.exit()
            srccur, srcconn = dbconnect(dbname, load_ext=True)
        except KeyError:
            print('If database is True, a file must be provided.')
        try:
            imglist = kwargs['images']
            # Select sources from one image at a time
            for img in imglist:
                srccur.execute('''SELECT id, bmaj, obs_ra, obs_dec, imsize, 
                    pixel_scale FROM Image WHERE name = ?''', (img, ))
                image = srccur.fetchone()
                srccur.execute('SELECT * FROM Source WHERE image_id = ?',
                               (image['id'], ))
                sources = srccur.fetchall()
                matches, non_matches, cat_matches = match(image, sources,
                                                          catalogs)
            srccur.close()
        except KeyError:
            print('If database is True, a list of images must be provided.')
    else:
        try:
            image, sources = kwargs['objects']
            # Make sure we're dealing with dictionaries for consistency
            # across data input types
            try:
                image['bmaj']
                dimage = image
            except TypeError:
                dimage = image.__dict__
            dsources = []
            for src in sources:
                try:
                    src['ra']
                    dsources.append(src)
                except TypeError:
                    dsources.append(src.__dict__)
            matches, non_matches, cat_matches = match(dimage, dsources,
                                                      catalogs)
        except KeyError:
            print('If database is False, an ImageTable object and list of '
                  'DetectedSource objects must be provided.')

    return matches, non_matches, cat_matches
