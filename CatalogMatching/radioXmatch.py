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
import catalogio
import crossmatch


def dbconnect(dbname):
    """Creates an sqlite3 cursor object and specifies
    row_factory so that fetch commands return the rows
    as dictionaries."""
    conn = sqlite3.connect(dbname)
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
    print('Adding catalog cross-match results to database Source table.\n')
    cur, conn = dbconnect(database)
    cur.execute('SELECT id FROM Image WHERE name = ?', (imname, ))
    img_id = cur.fetchone()[0]

    for src in matches + non_matches:
        cur.execute('''UPDATE Source SET catalog_id = ?, match_id = ?,
            min_deRuiter = ? WHERE src_id = ? AND image_id = ?''',
                    (src['catalog_id'], src['match_id'], src['min_deRuiter'],
                     src['src_id'], img_id))
    conn.commit()
    cur.close()


def search_range(image):
    """Returns the minimum & maximum RA, Dec of an image (in degrees).
    Used to limit the number of catalog sources extracted for 
    cross-matching."""
    try:
        # image must be a dictionary
        imsize = float(image['imsize'][1:].split(',')[0])
        fov = (image['pixel_scale'] * imsize) / 3600.
        fov = Angle(((image['pixel_scale'] * imsize) / 3600.) * u.deg)
    except TypeError:
        # imsize and/or pixel_scale are None
        fov = Angle(10. * u.deg) # default to 10 degrees
    try:      
        ra_ang = Angle(image['obs_ra'] * u.deg)
        ramin = (ra_ang - fov).wrap_at(360 * u.deg).degree
        ramax = (ra_ang + fov).wrap_at(360 * u.deg).degree
    except TypeError:
        # obs_ra is None
        ramin = None
        ramax = None
    try:
        dec_ang = Angle(image['obs_dec'] * u.deg)
        # ignore +/- 90 deg limit - won't matter
        decmin = (dec_ang - fov).degree
        decmax = (dec_ang + fov).degree
    except TypeError:
        # obs_dec is None
        decmin = None
        decmax = None
            
    search_range = (ramin, ramax, decmin, decmax)
    return search_range


def catalog_extract(cursor, catalogs, limits):
    """Returns a dictionary containing catalogs and their sources
    which lie in the specified range of RA, Dec."""
    ramin, ramax, decmin, decmax = limits
    print('Extracting sources from sky survey catalogs with RA between '
          '{} and {} and Dec between {} and {}.\n'.format(ramin, ramax,
                                                          decmin, decmax))
    catdict = {}
    for catalog in catalogs:
        if np.all(limits) is not None: # use both RA & Dec limits
            if ramin > ramax: # passes through 0
                query = 'SELECT * FROM %s WHERE ra BETWEEN ? AND ? AND dec '\
                        'BETWEEN ? AND ? OR ra BETWEEN ? AND ? AND dec '\
                        'BETWEEN ? AND ?' % catalog
                cursor.execute(query, (ramin, 360., decmin, decmax,
                                       0., ramax, decmin, decmax))
            else:
                query = 'SELECT * FROM %s WHERE ra BETWEEN ? AND ? AND dec '\
                        'BETWEEN ? AND ?' % catalog
                cursor.execute(query, (ramin, ramax, decmin, decmax))
            catdict[catalog] = cursor.fetchall()
        elif ramin is not None: # use only RA limits
            if ramin > ramax:
                query = 'SELECT * FROM %s WHERE ra BETWEEN ? AND ? AND ra '\
                        'BETWEEN ? AND ?' % catalog
                cursor.execute(query, (ramin, 360., 0., ramax))
            else:
                query = 'SELECT * FROM %s WHERE ra BETWEEN ? AND ?' % catalog
                cursor.execute(query, (ramin, ramax))
            catdict[catalog] = cursor.fetchall()
        else: # use only Dec limits
            query = 'SELECT * FROM %s WHERE dec BETWEEN ? AND ?' % catalog
            cursor.execute(query, (decmin, decmax))
            catdict[catalog] = cursor.fetchall()

    return catdict


def match(image, sources, catalogs):              
    """Performs the steps necessary to cross-match a 
    provided list of sources from the specified image
    to a given list of catalogs. The list of catalog
    names specifies which tables of the SkyCatalogs
    database to look in and in which order. The image
    is used to define a search range when extracting
    sources from the SkyCatalog database tables. The
    sources are cross-matched using the deRuiter
    radius criterion. A list of matched sources,
    non-matched sources, and the matched catalog
    sources is returned."""
    # Connect the sky catalog database
    skycat = os.path.join(catalogio.catalogdir, 'SkyCatalogs.sqlite')
    catcur, catconn = dbconnect(skycat)

    # Extract catalog sources
    catdict = catalog_extract(catcur, catalogs, search_range(image))

    # Define beam major axis FWHM:
    bmaj = image['bmaj'] / 3600. # deg

    # Cross-match each source to the catalog sources
    print('Searching for catalog cross-matches...\n')
    matches = []
    non_matches = []
    cat_matches = []
    for src in sources:
        idx = 0
        match, catsrc = crossmatch.deRuitermatch(
            src, catdict[catalogs[idx]], bmaj)
        while match != 1:
            idx += 1
            if idx < len(catalogs):
                match, catsrc = crossmatch.deRuitermatch(
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
            
    print('Matched {}/{} sources.\n'.format(len(matches), len(sources)))
    catcur.close()

    return matches, non_matches, cat_matches


def main(catalogs=['NVSS'], **kwargs):
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
    #      file='/data3/erichards/vlite/allvlite.sqlite')
    if kwargs['database'] is True:
        try:
            dbname = kwargs['file']
            if not os.path.exists(dbname):
                print('ERROR: File or path {} does not exist!'.format(dbname))
                sys.exit(0)
            else: pass
            srccur, srcconn = dbconnect(dbname)
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


if __name__ == '__main__':
         main()
