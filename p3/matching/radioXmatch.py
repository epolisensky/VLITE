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
import sqlite3
import numpy as np
from astropy.coordinates import Angle
import astropy.units as u
from skycatalog import catalogio
import matchfuncs


def dbconnect(dbname, load_ext=False):
    """Creates an sqlite3 cursor object and specifies
    row_factory so that fetch commands return the rows
    as dictionaries."""
    conn = sqlite3.connect(dbname)
    if load_ext:
        conn.enable_load_extension(True)
        conn.load_extension("./extension-functions.so")
        conn.enable_load_extension(False)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    return cur, conn


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


def cone_search(dbname, tables, center_coords, search_radius,
                same_res=False, **kwargs):
    cur, conn = dbconnect(dbname, load_ext=True)
    cra = center_coords[0]
    cdec = center_coords[1]
    r = np.radians(search_radius)
    d2r = np.pi / 180.
    print('\nExtracting sources from {} within {} degrees of {:.3f}, '
          '{:.3f}'.format(tables, int(round(search_radius)), cra, cdec))
    tabdict = {}
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

    if not same_res:
        for table in tables:
            query = '''SELECT * FROM %s WHERE 
                2 * ASIN(SQRT(SIN(((?-dec)/2)*?) * SIN(((?-dec)/2)*?) 
                + COS(?*?) * COS(dec*?) * SIN(((?-ra)/2)*?) 
                * SIN(((?-ra)/2)*?))) <= ?''' % table
            cur.execute(query, (cdec, d2r, cdec, d2r, cdec, d2r, d2r, cra, d2r,
                                cra, d2r, r))
            tabdict[table] = cur.fetchall()

    cur.close()

    return tabdict


#def match(image, sources, tabdict):


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
