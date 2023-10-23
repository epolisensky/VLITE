"""This module contains all the messy I/O functionality
for reading in the different sky survey catalogs.
Each catalog source is initialized as a CatalogSource
object to create a uniform set of properties for each
radio catalog's sources. When reading the catalogs for the
first time, sources with the uniform set of attributes
are written to text files called "[catalog_name]_psql.txt",
which is used to insert the sources into the PostgreSQL
database in an efficient way. 

Instructions for adding a new catalog:

1. Add catalog name in all lowercase and a new id number to 
the catalog dictionary. *DO NOT* re-order the id numbers 
unless you plan on re-writing all the "\*_psql.txt" files 
and the entire "radcat" schema. Otherwise, the id number of 
the catalog in the **radcat.catalogs** table will not match 
the 'catalog_id' of the sources in the catalog's table.

2. Write a function called ``read_[catalog name]`` for reading 
in the catalog and writing to the text file if it doesn't 
exist using the same format as all the others.

3. Add a line to radcatdb.py to call the 
``read_[catalog name]`` function.

Adapted from EP's iofuncs.py.

"""
import os,sys
import logging
import pandas as pd


# create logger
catio_logger = logging.getLogger('vdp.radiocatalogs.catalogio')


# Define path to the sky catalog data
catalogdir = '/home/vpipe/VLITE/VLITE/vdp/resources/RadioCatalogs'

# If adding a new catalog, don't forget to add it to this dictionary!
# Note: catalog name must be lowercase, cannot lead with a number, and
# cannot contain "."
catalog_dict = {'cosmos' : {'id' : 1}, 'first' : {'id' : 2},
                'gleam' : {'id' : 3}, 'gpsr1' : {'id' : 4},
                'gpsr5' : {'id' : 5}, 'lazio04' : {'id' : 6},
                'lofar_hba' : {'id' : 7}, 'lofar_lba' : {'id' : 8},
                'lotss' : {'id' : 9}, 'm31_glg04' : {'id' : 10},
                'nordgc' : {'id' : 11}, 'nrl_nvss' : {'id' : 12},
                'nvss' : {'id' : 13}, 'sevenc' : {'id' : 14},
                'sumss' : {'id' : 15}, 'tgss' : {'id' : 16},
                'txs' : {'id' : 17}, 'vlssr' : {'id' : 18},
                'wenss': {'id' : 19}, 'atlas': {'id' : 20},
                'lotssdr1' : {'id' : 21}, 'racs' : {'id' : 22},
                'psr' : {'id' : 23}, 'vlass1' : {'id' : 24},
                'vlass2' : {'id' : 25}}
                #'vlass2' : {'id' : 25}, 'vlass3' : {'id' : 26}}


def dms2deg(d, m, s):
    """Translates coordinates from deg:min:sec to decimal degrees.
    If computing RA, the result needs to be multiplied by 15.

    """
    dec = abs(float(d)) + abs(float(m))/60. + abs(float(s))/3600.
    d = str(d)
    if d.__contains__("-"):
        dec = -dec
    return float(dec)


def set_error(catalog, attrs):
    """Sets missing values for specified attributes to
    the median value of the distribution.

    Parameters
    ----------
    catalog : list
        List of CatalogSource objects with missing
        attribute values.
    attrs : list
        Attribute names which contain missing values to be
        set to the median of the distribution.
   
    Returns
    -------
    catalog : list
        CatalogSource objects with modified attributes
        that contain no missing values.
    """
    for attr in attrs:
        s = pd.Series([getattr(src, attr) for src in catalog])
        s.loc[s.notnull() == False] = s.median()
        for idx in range(len(catalog)):
            setattr(catalog[idx], attr, s[idx])
    return catalog


class CatalogSource(object):
    """Class for the radio catalog sources. Their
    catalog origin is identified by the 'catalog_id'
    attribute.

    Attributes
    ----------
    id : int
        Uniquely identifies the source in its catalog.
    name : str
        Name given to the source in some catalogs.
    ra : float
        Source right ascension (degrees).
    e_ra : float
        Error on the right ascension (degrees).
    dec : float
        Source declination (degrees).
    e_dec : float
        Error on the declination (degrees).
    total_flux : float
        Total integrated flux (mJy).
    e_total_flux : float
        Error on the total flux (mJy).
    peak_flux : float
        Peak flux density per beam (mJy/beam).
    e_peak_flux : float
        Error on the peak flux (mJy/beam).
    maj : float
        FWHM of the source major axis (arcsec).
    e_maj : float
        Error on the major axis size (arcsec).
    min : float
        FWHM of the source minor axis (arcsec).
    e_min : float
        Error on the minor axis size (arcsec).
    pa : float
        Position angle of the source major axis (degrees).
    e_pa : float
        Error on the position angle (degrees).
    rms : float
        Local or image noise (mJy/beam).
    field : str
        Name of the field where the source is located
        provided by some catalogs.
    catalog_id : int
        Uniquely identifies the radio catalog from which this
        source originates.
    pt_like : boolean
        True if source is point-like. Null if unknown
    assoc_id : int
        The id of the closest VLITE source in the
        **assoc_source** table.
    sep : float
        Angular separation between this catalog source and the
        closest VLITE source (arcsec).
    duplicate_flag : float
        For VLASS catalogs. Select < 2 for high quality sources
    quality_flag : float
        For VLASS catalogs. Select =0 or =4 for high quality sources
    p_sidelobe : float
        For VLASS catalogs. Probability source is a sidelobe
    """
    def __init__(self):
        self.id = None
        self.name = None
        self.ra = None # deg
        self.e_ra = None # deg
        self.dec = None # deg
        self.e_dec = None # deg
        self.total_flux = None # mJy
        self.e_total_flux = None # mJy
        self.peak_flux = None # mJy/beam
        self.e_peak_flux = None # mJy/beam
        self.maj = None # arcsec
        self.e_maj = None # arcsec
        self.min = None # arcsec
        self.e_min = None # arcsec
        self.pa = None # deg
        self.e_pa = None # deg
        self.rms = None # mJy/beam
        self.field = None
        self.catalog_id = None
        self.pt_like = None
        self.assoc_id = None
        self.sep = None
        self.duplicate_flag = None
        self.quality_flag = None
        self.p_sidelobe = None

def read_tgss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the TGSS survey catalog and writes them into a file in the 
    same directory called tgss_psql.txt if the file does
    not already exist.

    Telescope/frequency: GMRT 150 MHz

    Spatial resolution: 25'' 

    """
    catalog_dict['tgss']['telescope'] = 'GMRT'
    catalog_dict['tgss']['frequency'] = 150.
    catalog_dict['tgss']['resolution'] = 25.
    catalog_dict['tgss']['reference'] = 'Intema et al. (2017)'
    
    psqlf = os.path.join(catalogdir, 'tgss_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir, 'TGSSADR1_7sigma_catalog.tsv')
    f = open(fname, 'r')
    b = f.readlines()
    f.close()    
    n = len(b)
    cnt = 1
    sources = []
    for i in range(1, n):
        sources.append(CatalogSource())
        d = b[i].split() # Data Line
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = d[0]+'_'+d[1]
        sources[-1].ra = float(d[2]) # deg
        sources[-1].e_ra = float(d[3])/3600.0 # deg
        sources[-1].dec = float(d[4]) # deg
        sources[-1].e_dec = float(d[5])/3600.0 # deg
        sources[-1].total_flux = float(d[6]) # mJy
        sources[-1].e_total_flux = float(d[7]) # mJy
        sources[-1].peak_flux = float(d[8]) # mJy/beam
        sources[-1].e_peak_flux = float(d[9]) # mJy/beam
        sources[-1].maj = float(d[10]) # arcsec
        sources[-1].e_maj = float(d[11]) # arcsec
        sources[-1].min = float(d[12]) # arcsec
        sources[-1].e_min = float(d[13]) # arcsec
        sources[-1].pa = float(d[14]) # deg
        sources[-1].e_pa = float(d[15]) # deg
        sources[-1].rms = float(d[16]) # mJy/beam
        sources[-1].Code = d[17]
        sources[-1].field = d[18]
        sources[-1].catalog_id = catalog_dict['tgss']['id']
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} TGSS sources to tgss_psql.txt'.format(
        len(sources)))
    return sources


def read_first(return_sources=False):
    """Generates a list of CatalogSource objects from
    the FIRST survey catalog and writes them into a file in the 
    same directory called first_psql.txt if the file does
    not already exist.

    Telescope/frequency: VLA 1.4 GHz

    Spatial resolution: 5'' 

    """
    catalog_dict['first']['telescope'] = 'VLA'
    catalog_dict['first']['frequency'] = 1400.
    catalog_dict['first']['resolution'] = 5.
    catalog_dict['first']['reference'] = 'White et al. (1997)'
    
    psqlf = os.path.join(catalogdir, 'first_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    sources = []
    fname = os.path.join(catalogdir, 'FIRST_catalog_14dec17.txt')
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        if line[0] != '#':
            line = line.split()
            if float(line[6]) < (0.15):
                sources.append(CatalogSource())
                sources[-1].id = cnt
                sources[-1].name = 'FIRST_'+str(cnt)
                cnt += 1
                sources[-1].ra = 15. * dms2deg(line[0], line[1], line[2]) # deg
                sources[-1].dec = dms2deg(line[3], line[4], line[5]) # deg
                sources[-1].prob = float(line[6])
                sources[-1].peak_flux = float(line[7]) # mJy/beam
                sources[-1].total_flux = float(line[8]) # mJy
                sources[-1].rms = float(line[9]) # mJy/beam
                sources[-1].maj = float(line[10]) # arcsec
                sources[-1].min = float(line[11]) # arcsec
                sources[-1].pa = float(line[12])
                # HACK: FIRST positional uncertainty calculations given
                # for FITTED bmaj & bmin. 
                #  I'll assign these to e_ra & e_dec for lack of better method
                snr = (sources[-1].peak_flux-0.25) / sources[-1].rms
                dfMaj = float(line[13]) * ((1.0/snr) + 0.05) / 3600.0 # deg
                dfMin = float(line[14]) * ((1.0/snr) + 0.05) / 3600.0 # deg
                # fpa = float(line[15]) # deg
                sources[-1].e_ra = dfMin # deg
                sources[-1].e_dec = dfMaj # deg
                sources[-1].catalog_id = catalog_dict['first']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} FIRST sources to first_psql.txt'.format(
        len(sources)))
    return sources


def read_sumss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the SUMSS survey catalog and writes them into a file in the 
    same directory called sumss_psql.txt if the file does
    not already exist.

    Telescope/frequency: MOST 843 MHz

    Spatial resolution: 45'' 

    """
    catalog_dict['sumss']['telescope'] = 'MOST'
    catalog_dict['sumss']['frequency'] = 843.
    catalog_dict['sumss']['resolution'] = 45.
    catalog_dict['sumss']['reference'] = 'Mauch et al. (2003)'
    
    psqlf = os.path.join(catalogdir, 'sumss_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    sources = []
    fname = os.path.join(catalogdir, 'SUMSS.txt')
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        if line[0] != ';':
            sources.append(CatalogSource())
            sources[-1].id = cnt
            sources[-1].name = 'SUMSS_'+str(cnt)
            cnt += 1
            line = line.split()
            sources[-1].ra = 15. * dms2deg(line[0], line[1], line[2]) # deg
            sources[-1].dec = dms2deg(line[3], line[4], line[5]) # deg
            sources[-1].e_ra = float(line[6])/3600.0 # deg
            sources[-1].e_dec = float(line[7])/3600.0 # deg
            sources[-1].peak_flux = float(line[8]) # mJy/bm
            sources[-1].e_peak_flux = float(line[9])
            sources[-1].total_flux = float(line[10]) # mJy
            sources[-1].e_total_flux = float(line[11])
            sources[-1].maj = float(line[12]) # arcsec
            sources[-1].min = float(line[13]) # arcsec
            sources[-1].pa = float(line[14])
            sources[-1].catalog_id = catalog_dict['sumss']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} SUMSS sources to sumss_psql.txt'.format(
        len(sources)))
    return sources


def read_wenss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the WENSS survey catalog and writes them into a file in the 
    same directory called wenss_psql.txt if the file does
    not already exist.

    Telescope/frequency: WSRT 325 MHz

    Spatial resolution: 54'' 
    
    """
    catalog_dict['wenss']['telescope'] = 'WSRT'
    catalog_dict['wenss']['frequency'] = 325.
    catalog_dict['wenss']['resolution'] = 54.
    catalog_dict['wenss']['reference'] = 'Rengelink et al. (1997)'
    
    psqlf = os.path.join(catalogdir, 'wenss_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    sources = []
    fname = os.path.join(catalogdir, 'WENS.COMPLETE')
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        if line[0] != ';':
            sources.append(CatalogSource())
            line = line.split()
            sources[-1].id = cnt
            sources[-1].name = 'WENSS_'+str(cnt)
            cnt += 1
            sources[-1].ra = float(line[0]) # deg
            sources[-1].dec = float(line[1]) # deg
            sources[-1].total_flux = float(line[2]) # mJy
            sources[-1].e_total_flux = float(line[3]) # mJy
            sources[-1].peak_flux = float(line[4]) # mJy/beam
            sources[-1].e_peak_flux = float(line[5]) # mJy/beam
            sources[-1].rms = float(line[6]) # mJy/beam
            sources[-1].e_ra = float(line[7])/3600.0 # deg
            sources[-1].e_dec = float(line[8])/3600.0 # deg
            sources[-1].maj = float(line[9]) # arcsec
            sources[-1].min = float(line[10]) # arcsec
            sources[-1].pa = float(line[11])
            sources[-1].pt_like = line[12]
            sources[-1].catalog_id = catalog_dict['wenss']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} WENSS sources to wenss_psql.txt'.format(
        len(sources)))
    return sources


def read_nvss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the NVSS survey catalog and writes them into a file in the 
    same directory called nvss_psql.txt if the file does
    not already exist.

    Telescope/frequency: VLA 1.4 GHz

    Spatial resolution: 45'' 

    """
    catalog_dict['nvss']['telescope'] = 'VLA'
    catalog_dict['nvss']['frequency'] = 1400.
    catalog_dict['nvss']['resolution'] = 45.
    catalog_dict['nvss']['reference'] = 'Condon et al. (1998)'
    
    psqlf = os.path.join(catalogdir, 'nvss_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    sources = []
    fname = os.path.join(catalogdir, 'NVSSCatalog.text')
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line:
            break
        beginning = line[0:1]
        later = line[7:8]
        if beginning.isdigit() == True:
            sources.append(CatalogSource())
            sources[-1].id = cnt
            sources[-1].name = 'NVSS_'+str(cnt)
            cnt += 1
            if line[6] == ' ': # ra
                sources[-1].ra = 15. * dms2deg(line[0:2], line[3:5],
                                               line[7:11])
            else:
                sources[-1].ra = 15. * dms2deg(line[0:2], line[3:5],
                                               line[6:11])
            if line[19] == ' ': # dec
                sources[-1].dec = dms2deg(line[12:15], line[16:18],
                                          line[20:23])
            else:
                sources[-1].dec = dms2deg(line[12:15], line[16:18],
                                          line[19:23])
            if line[25] == ' ' and line[26] == ' ': # flux
                sources[-1].total_flux = float(line[27:30]) # mJy
            elif line[25] == ' ' and line[26] != ' ':
                sources[-1].total_flux = float(line[26:30])
            else:
                sources[-1].total_flux = float(line[25:30])
            if line[32] == ' ': #maj
                sources[-1].maj  = float(line[33:36]) # arcsec
            elif line[31] == ' ' and line[32] != ' ':
                sources[-1].maj = float(line[32:36])
            elif line[31]=='<':
                sources[-1].maj = float(line[32:36])
            else:
                sources[-1].maj = float(line[31:36])
            if line[37] == '<': # min
                sources[-1].min = float(line[38:42]) # arcsec
            elif line[37] == ' ' and line[38] != ' ':
                sources[-1].min = float(line[38:42])
            else:
                sources[-1].min = float(line[39:42])
            if line[43] == '-': # pa
                sources[-1].pa = float(line[43:48])
            elif line[43] == ' ' and line[44] == ' ' and line[45] != ' ':
                sources[-1].pa = float(line[45:48])
            elif line[43] == ' ' and line[44] != ' ':
                sources[-1].pa = float(line[44:48])
            else:
                sources[-1].pa = 0.0
            sources[-1].field = line[65:73]
        if beginning.isdigit() == False and line[0:7] == '       ':
            sources[-1].e_ra = float(15*dms2deg(0,0,line[7:11])) # error ra
            if line[19] == ' ': # error dec
                sources[-1].e_dec = float(dms2deg(0,0,line[20:23]))
            else:
                sources[-1].e_dec = float(dms2deg(0,0,line[19:23]))
            sources[-1].e_total_flux = float(line[27:30]) # error flux (mJy)
            sources[-1].catalog_id = catalog_dict['nvss']['id']
    fread.close()
    sources = set_error(sources, ['e_ra', 'e_dec'])
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %s %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} NVSS sources to nvss_psql.txt'.format(
        len(sources)))
    return sources


def read_gleam(return_sources=False):
    """Generates a list of CatalogSource objects from
    the GLEAM survey catalog and writes them into a file in the 
    same directory called gleam_psql.txt if the file does
    not already exist.

    Telescope/frequency: MWA 74-231 MHz

    Spatial resolution: ~100'' 

    """
    catalog_dict['gleam']['telescope'] = 'MWA'
    catalog_dict['gleam']['frequency'] = 150.
    catalog_dict['gleam']['resolution'] = 100.
    catalog_dict['gleam']['reference'] = 'Hurley-Walker et al. (2017)'
    
    psqlf = os.path.join(catalogdir, 'gleam_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir, 'gleamegc.dat')
    f = open(fname, 'r')
    b = f.readlines()
    f.close()    
    n = len(b)
    #print 'reading %s' % fname
    #print ' read %d lines' % n
    cnt = 1
    sources = []
    for i in range(0, n):
        sources.append(CatalogSource())
        d = b[i].split() # Data Line
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = d[0]+'_'+d[1]
        sources[-1].ra = float(d[10]) # deg
        if d[11].replace('.','',1).isdigit():
            sources[-1].e_ra = float(d[11]) # deg
        sources[-1].dec = float(d[12]) # deg
        if d[13].replace('.','',1).isdigit():
            sources[-1].e_dec = float(d[13]) # deg
        sources[-1].total_flux = float(d[16])*1000.0 # mJy
        sources[-1].e_total_flux = float(d[17])*1000.0 # mJy
        sources[-1].peak_flux = float(d[14])*1000.0 # mJy/beam
        sources[-1].e_peak_flux = float(d[15])*1000.0 # mJy/beam
        sources[-1].maj = float(d[18]) # arcsec
        sources[-1].e_maj = float(d[19]) # arcsec
        sources[-1].min = float(d[20]) # arcsec
        sources[-1].e_min = float(d[21]) # arcsec
        sources[-1].pa = float(d[22]) # deg
        sources[-1].e_pa = float(d[23]) # deg
        # sources[-1].rms = float(d[16]) # mJy/beam
        sources[-1].catalog_id = catalog_dict['gleam']['id']
    sources = set_error(sources, ['e_ra', 'e_dec'])
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} GLEAM sources to gleam_psql.txt'.format(
        len(sources)))
    return sources


def read_cosmos(return_sources=False):
    """Generates a list of CatalogSource objects from
    the COSMOS Legacy survey catalog and writes them 
    into a file in the same directory called cosmos_psql.txt
    if the file does not already exist. All non-header
    lines start with 'C'.

    Telescope/frequency: VLA 320 MHz

    Spatial resolution: ~6''

    """
    catalog_dict['cosmos']['telescope'] = 'VLA'
    catalog_dict['cosmos']['frequency'] = 320.
    catalog_dict['cosmos']['resolution'] = 6.
    catalog_dict['cosmos']['reference'] = 'Smolcic et al. (2014)'
    
    psqlf = os.path.join(catalogdir, 'cosmos_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir,
                         'vla-cosmos_327_sources_published_version.tbl')
    fin = open(fname, 'r')
    cnt = 1
    sources = []
    # Read file:
    while 1:
        line2 = fin.readline()
        # If end of file, break out of loop:
        if not line2: break
        # Check if data line:
        if line2[0] == 'C':
            sources.append(CatalogSource())
            line = line2.split()
            sources[-1].id = cnt
            cnt += 1
            sources[-1].name = line[0] # name in 90cm catalog
            sources[-1].ra = float(line[1]) # deg
            sources[-1].dec = float(line[2]) # deg
            sources[-1].e_ra = float(line[9])/3600.0 # deg
            sources[-1].e_dec = float(line[10])/3600.0 # deg
            sources[-1].peak_flux = float(line[11]) # mJy/beam
            sources[-1].e_peak_flux = float(line[12]) # mJy/beam
            sources[-1].total_flux = float(line[13]) # mJy
            sources[-1].e_total_flux = float(line[14]) # mJy
            sources[-1].rms = float(line[15]) # mJy/beam       
            sources[-1].maj = float(line[16]) # arcsec
            sources[-1].min = float(line[17]) # arcsec
            sources[-1].pa = float(line[18]) # deg
            # 0 = unresolved, 1 =resolved
            sources[-1].flagresolved = int(line[19])
            if int(line[19])==0: sources[-1].pt_like = 'True'
            else: sources[-1].pt_like = 'False'
            # 0 = single component, 1 = multi-component
            sources[-1].flagmulticom = int(line[20])
            sources[-1].name20 = line[21] # name in 20cm (1.4 GHz) catalog
            sources[-1].Peak20 = float(line[22]) # mJy/beam, 20cm
            sources[-1].Total20 = float(line[23]) # mJy, 20cm
            # separation between 90cm and 20cm positions 
            if line[24][0]!='*':
                sources[-1].sep = float(line[24])/3600.0 # deg
            if sources[-1].maj < 1.: sources[-1].maj = 1.
            if sources[-1].min < 1.: sources[-1].min = 1.
            sources[-1].catalog_id = catalog_dict['cosmos']['id']
    fin.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} COSMOS Legacy P-band sources to '
                      'cosmos_psql.txt'.format(len(sources)))
    return sources


def read_vlssr(return_sources=False):
    """Generates a list of CatalogSource objects from
    the VLA Low-frequency Sky Survey Redux catalog and
    writes them into a file in the same directory called
    vlssr_psql.txt if the file does not already exist.

    Telescope/frequency: VLA 74 MHz

    Spatial resolution: ~75''

    """
    catalog_dict['vlssr']['telescope'] = 'VLA'
    catalog_dict['vlssr']['frequency'] = 74.
    catalog_dict['vlssr']['resolution'] = 75.
    catalog_dict['vlssr']['reference'] = 'Lane et al. (2014)'
    
    psqlf = os.path.join(catalogdir, 'vlssr_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir, 'FullVLSSCatalog.clean')
    fread = open(fname, 'r')
    sources = []
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        sources.append(CatalogSource())
        sources[-1].id = cnt
        sources[-1].name = 'VLSSr_'+str(cnt)
        cnt += 1
        sources[-1].ra = 15. * dms2deg(line[0:2], line[3:5], line[6:11]) # deg
        sources[-1].dec = dms2deg(line[12:15], line[16:18], line[19:23]) # deg
        sources[-1].total_flux = float(line[30:36]) * 1000.0 # mJy
        sources[-1].maj = float(line[38:42]) # arcsec
        sources[-1].min = float(line[44:48]) # arcsec
        sources[-1].pa = float(line[49:54]) # deg
        # Read errors line
        line = fread.readline()
        if not line:
            catio_logger.error('ERROR: Cannot read error line in VLSSr catalog')
            break
        sources[-1].e_ra = 15. * dms2deg(0, 0, line[6:11]) # deg
        sources[-1].e_dec = dms2deg(0, 0, line[19:23]) # deg
        sources[-1].e_total_flux = float(line[30:36]) * 1000.0 # mJy
        sources[-1].catalog_id = catalog_dict['vlssr']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} VLSSr sources to vlssr_psql.txt'.format(
        len(sources)))
    return sources


def read_txs(return_sources=False):
    """Generates a list of CatalogSource objects from
    the TXS survey catalog and writes them into a file in the 
    same directory called txs_psql.txt if the file does
    not already exist.

    Telescope/frequency: Texas Interferometer 365 MHz

    Spatial resolution: 100"? 

    """
    catalog_dict['txs']['telescope'] = 'Texas Interferometer'
    catalog_dict['txs']['frequency'] = 365.
    catalog_dict['txs']['resolution'] = 100.
    catalog_dict['txs']['reference'] = 'Douglas et al. (1996)'
    
    psqlf = os.path.join(catalogdir, 'txs_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir, 'TXS_J2000.txt')
    fread = open(fname, 'r')
    # First line is header
    line = fread.readline()
    sources = []
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        sources.append(CatalogSource())
        line = line.split()
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        sources[-1].ra = float(line[1]) # deg
        sources[-1].dec = float(line[2]) # deg
        sources[-1].e_ra = float(line[3]) # deg
        sources[-1].e_dec = float(line[4]) # deg
        sources[-1].total_flux = float(line[5]) # mJy
        sources[-1].e_total_flux = float(line[6]) # mJy
        sources[-1].catalog_id = catalog_dict['txs']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} TXS sources to txs_psql.txt'.format(
        len(sources)))
    return sources


def read_sevenc(return_sources=False):
    """Generates a list of CatalogSource objects from
    the 7C survey catalog and writes them into a file in the 
    same directory called sevenc_psql.txt if the file does
    not already exist.

    Telescope/frequency: Cambridge Low Frequency Synthesis Telescope 151 MHz

    Spatial resolution: ~70''

    """
    catalog_dict['sevenc']['telescope'] = 'CLFST'
    catalog_dict['sevenc']['frequency'] = 151.
    catalog_dict['sevenc']['resolution'] = 70.
    catalog_dict['sevenc']['reference'] = 'Riley at al. (1999)'
    
    psqlf = os.path.join(catalogdir, 'sevenc_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir, '7C_new.txt')
    fread = open(fname, 'r')
    sources = []
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        beg = line[0]
        if beg.isdigit():
            sources.append(CatalogSource())
            sources[-1].id = cnt
            sources[-1].name = '7C_'+str(cnt)
            cnt += 1
            line = line.split()
            sources[-1].ra = float(line[0])
            sources[-1].e_ra = float(line[1])/3600.0
            sources[-1].dec = float(line[2])
            sources[-1].e_dec = float(line[3])/3600.0
            sources[-1].peak_flux = float(line[4])
            sources[-1].e_peak_flux = float(line[5])
            sources[-1].total_flux = float(line[6])
            sources[-1].e_total_flux = float(line[7])
            sources[-1].snr = float(line[8])
            sources[-1].catalog_id = catalog_dict['sevenc']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} 7C sources to sevenc_psql.txt'.format(
        len(sources)))
    return sources


def read_gpsr5(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Galactic Plan 5-GHz VLA Survey catalog and writes
    them into a file in the same directory called gpsr5_psql.txt
    if the file does not already exist.

    Telescope/frequency: VLA 5 GHz

    Spatial resolution: ~4''

    """
    catalog_dict['gpsr5']['telescope'] = 'VLA'
    catalog_dict['gpsr5']['frequency'] = 5000.
    catalog_dict['gpsr5']['resolution'] = 4.
    catalog_dict['gpsr5']['reference'] = 'Becker et al. (1994)'
    
    psqlf = os.path.join(catalogdir, 'gpsr5_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir, 'GPRS_5GHz.txt')
    fread = open(fname, 'r')
    sources = []
    # Skip header line
    line = fread.readline()
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        sources[-1].ra = float(line[1]) # deg
        sources[-1].dec = float(line[2]) # deg
        sources[-1].e_ra = float(line[3]) # deg
        sources[-1].e_dec = float(line[4]) # deg
        sources[-1].peak_flux = float(line[5]) # mJy/beam
        sources[-1].total_flux = float(line[6]) # mJy
        sources[-1].size = float(line[7]) # diameter, arcsec
        if sources[-1].size < 1.:
            sources[-1].maj = 0.5 # arcsec
            sources[-1].min = 0.5 # arcsec 
        else:
            sources[-1].maj = sources[-1].size * 0.5 # arcsec
            sources[-1].min = sources[-1].size * 0.5 # arcsec
        sources[-1].pa = 0.0 # deg
        sources[-1].catalog_id = catalog_dict['gpsr5']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} GPSR 5 GHz sources to gpsr5_psql.txt'.
                      format(len(sources)))
    return sources


def read_gpsr1(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Galactic Plan 1.4-GHz VLA Survey catalog and writes
    them into a file in the same directory called gpsr1_psql.txt
    if the file does not already exist.

    Telescope/frequency: VLA 1.4 GHz

    Spatial resolution: ~5''

    """
    catalog_dict['gpsr1']['telescope'] = 'VLA'
    catalog_dict['gpsr1']['frequency'] = 1400.
    catalog_dict['gpsr1']['resolution'] = 5.
    catalog_dict['gpsr1']['reference'] = 'Zoonematkermani et al. (1990)'
    
    psqlf = os.path.join(catalogdir, 'gpsr1_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    fname = os.path.join(catalogdir, 'GPRS_1.4GHz.txt')
    fread = open(fname, 'r')
    sources = []
    # Skip header line
    line = fread.readline()
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        sources[-1].ra = float(line[1]) # deg
        sources[-1].dec = float(line[2]) # deg
        sources[-1].e_ra = float(line[3]) # deg
        sources[-1].e_dec = float(line[4]) # deg
        sources[-1].peak_flux = float(line[5]) # mJy/beam
        sources[-1].total_flux = float(line[6]) # mJy
        sources[-1].size = float(line[7]) # diameter, arcsec
        if sources[-1].size < 1.:
            sources[-1].maj = 0.5 # arcsec
            sources[-1].min = 0.5 # arcsec 
        else:
            sources[-1].maj = sources[-1].size * 0.5 # arcsec
            sources[-1].min = sources[-1].size * 0.5 # arcsec
        sources[-1].pa = 0.0 # deg
        sources[-1].catalog_id = catalog_dict['gpsr1']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} GPSR 1.4 GHz sources to gpsr1_psql.txt'.
                      format(len(sources)))
    return sources


def read_nordgc(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Nord et al. Galactic Center survey catalog and writes
    them into a file in the same directory called nordgc_psql.txt
    if the file does not already exist.

    NOTE: This is actually the Hyman-updated version of the
    Nord et al. (2004) catalog.

    Telescope/frequency: VLA 330 MHz

    Spatial resolution: ~10''

    """
    catalog_dict['nordgc']['telescope'] = 'VLA'
    catalog_dict['nordgc']['frequency'] = 330.
    catalog_dict['nordgc']['resolution'] = 10.
    catalog_dict['nordgc']['reference'] = 'Nord et al. (2004)'
    
    psqlf = os.path.join(catalogdir, 'nordgc_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    #fname = os.path.join(catalogdir, 'NordGC_330MHz.txt')
    fname = os.path.join(catalogdir, 'NADAETC.startable')
    fread = open(fname, 'r')
    sources = []
    # Skip header line
    #line = fread.readline()
    # Skip 3 header lines
    line = fread.readline()
    line = fread.readline()
    line = fread.readline()
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        #sources[-1].name = line[0]
        sources[-1].name = 'HymanNordGC_%d' % cnt
        cnt += 1
        #sources[-1].ra = float(line[1]) # deg       
        sources[-1].ra = 15. * dms2deg(line[0], line[1], line[2]) # deg
        #sources[-1].dec = float(line[2]) # deg
        sources[-1].dec = dms2deg(line[3], line[4], line[5]) # deg
        #sources[-1].e_ra = float(line[3]) # deg
        sources[-1].e_ra = 5./3600. # approx 5 arcsec
        #sources[-1].e_dec = float(line[4]) # deg
        sources[-1].e_dec = 5./3600.
        #sources[-1].peak_flux = float(line[5]) # mJy/bm
        #sources[-1].rms = float(line[6]) # mJy/bm
        #sources[-1].total_flux = float(line[7]) # mJy
        #sources[-1].size = float(line[8]) # diameter, arcsec
        #if (sources[-1].size < 1.):
        #    sources[-1].maj = 0.5 # arcsec
        #    sources[-1].min = 0.5 # arcsec 
        #else:
        #    sources[-1].maj = sources[-1].size*0.5 # arcsec
        #    sources[-1].min = sources[-1].size*0.5 # arcsec
        #sources[-1].pa = 0.0 # deg
        # Assuming beam size for src sizes
        sources[-1].maj = 0.5 * 13. # arcsec
        sources[-1].min = 0.5 * 4.5 # arcsec
        sources[-1].pa = 0. # deg
        sources[-1].catalog_id = catalog_dict['nordgc']['id']
    fread.close()
    #print 'Read %d Nord GC 330 MHz catalog sources' % len(sources)
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} Hyman-updated Nord GC 330 MHz catalog '
                      'sources to nordgc_psql.txt'.format(len(sources)))
    return sources


def read_lazio04(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Nord et al. Galactic Center survey catalog and writes
    them into a file in the same directory called lazio04_psql.txt
    if the file does not already exist.

    Telescope/frequency: VLA 330 MHz

    Spatial resolution: 6?''

    """
    catalog_dict['lazio04']['telescope'] = 'VLA'
    catalog_dict['lazio04']['frequency'] = 330.
    catalog_dict['lazio04']['resolution'] = 6.
    catalog_dict['lazio04']['reference'] = 'Lazio (2004)'
    
    psqlf = os.path.join(catalogdir, 'lazio04_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir, 'Lazio2004.txt')
    fread = open(fname, 'r')
    sources = []
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        h = int(line[1])
        m = int(line[2])
        s = float(line[3])
        sources[-1].ra = 15. * dms2deg(h, m, s) # deg
        d = int(line[4])
        m = int(line[5])
        s = float(line[6])
        sources[-1].dec = dms2deg(d, m, s) # deg
        # uncertainties not given, assuming 1 arcsec
        sources[-1].e_ra = 1./3600. # deg
        sources[-1].e_dec = 1./3600. # deg
        sources[-1].maj = 0.5 # arcsec
        sources[-1].min = 0.5 # arcsec 
        sources[-1].pa = 0. # deg
        sources[-1].catalog_id = catalog_dict['lazio04']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} Lazio 2004 catalog sources to '
                      'lazio04_psql.txt'.format(len(sources)))
    return sources


def read_m31_glg04(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Nord et al. Galactic Center survey catalog and writes
    them into a file in the same directory called m31_glg04_psql.txt
    if the file does not already exist.

    Telescope/frequency: VLA 325 MHz

    Spatial resolution: 6''

    """
    catalog_dict['m31_glg04']['telescope'] = 'VLA'
    catalog_dict['m31_glg04']['frequency'] = 325.
    catalog_dict['m31_glg04']['resolution'] = 6.
    catalog_dict['m31_glg04']['reference'] = 'Gelfand, Lazio, & Gaensler (2004)'
    
    psqlf = os.path.join(catalogdir, 'm31_glg04_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir, 'GLG2004_M31_TABLE3.dat')
    sources = []
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        sources[-1].morph = line[1]
        h = int(line[2])
        m = int(line[3])
        s = float(line[4])
        sources[-1].ra = 15. * dms2deg(h, m, s) # deg
        d = int(line[6])
        m = int(line[7])
        s = float(line[8])
        sources[-1].dec = dms2deg(d, m, s) # deg
        sources[-1].e_ra = float(line[5])/3600.0 # deg
        sources[-1].e_dec = float(line[9])/3600.0 # deg
        sources[-1].maj = float(line[12]) # arcsec
        sources[-1].e_maj = float(line[13]) # arcsec
        sources[-1].min = float(line[14]) # arcsec
        sources[-1].e_min = float(line[15]) # arcsec
        sources[-1].pa = float(line[16]) # deg
        sources[-1].e_pa = float(line[17]) # deg
        sources[-1].peak_flux = float(line[18]) # mJy/bm
        sources[-1].e_peak_flux  = float(line[19]) # mJy/bm
        sources[-1].total_flux = float(line[20]) # mJy
        sources[-1].e_total_flux = float(line[21]) # mJy
        sources[-1].rms = float(line[22]) # mJy/bm
        sources[-1].catalog_id = catalog_dict['m31_glg04']['id']
        #if sources[-1].maj < 3e-4: 
        #    print '  %s %e' % (sources[-1].id, sources[-1].maj)
        #if sources[-1].min < 3e-4: 
        #    print '  %s %e' % (sources[-1].id, sources[-1].min)
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} GLG2004 M31 catalog sources to '
                      'm31_glg04_psql.txt'.format(len(sources)))
    return sources


def read_lotss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Nord et al. Galactic Center survey catalog and writes
    them into a file in the same directory called lotss_psql.txt
    if the file does not already exist.

    Telescope/frequency: LOFAR 150 MHz

    Spatial resolution: 25''

    """
    catalog_dict['lotss']['telescope'] = 'LOFAR'
    catalog_dict['lotss']['frequency'] = 150
    catalog_dict['lotss']['resolution'] = 25.
    catalog_dict['lotss']['reference'] = 'Shimwell et al. (2017)'
    
    psqlf = os.path.join(catalogdir, 'lotss_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir, 'lotss.dat')
    sources = []
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        sources[-1].ra = float(line[1]) # deg
        sources[-1].e_ra = float(line[3])/3600. # deg (total 1sigma error)
        sources[-1].dec = float(line[4]) # deg
        sources[-1].e_dec = float(line[6])/3600. # deg (total 1 sigma error)
        sources[-1].peak_flux = float(line[7]) # mJy/bm
        sources[-1].e_peak_flux = float(line[9]) # mJy/bm
        sources[-1].total_flux = float(line[10]) # mJy
        sources[-1].e_total_flux = float(line[12]) # mJy
        sources[-1].rflag = line[13]
        sources[-1].rms = float(line[14]) # mJy/bm
        sources[-1].code = line[15]
        sources[-1].field = line[16]
        # these are not given in the catalog, setting to half resolution:
        sources[-1].maj = 12.5 # arcsec
        sources[-1].min = 12.5 # arcsec
        sources[-1].pa = 0. # deg
        sources[-1].catalog_id = catalog_dict['lotss']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} LoTSS sources to lotss_psql.txt'.format(
        len(sources)))
    return sources


def read_lofar_lba(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Nord et al. Galactic Center survey catalog and writes
    them into a file in the same directory called lofar_lba_psql.txt
    if the file does not already exist.

    Telescope/frequency: LOFAR 34, 46, & 62 MHz

    Spatial resolution: 30-56''

    """
    catalog_dict['lofar_lba']['telescope'] = 'LOFAR'
    catalog_dict['lofar_lba']['frequency'] = 62
    catalog_dict['lofar_lba']['resolution'] = 40.
    catalog_dict['lofar_lba']['reference'] = 'van Weeren et al. (2014)'
    
    psqlf = os.path.join(catalogdir, 'lofar_lba_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir, 'LOFAR_LBA.dat')
    sources = []
    fread = open(fname, 'r')
    data_lines = fread.readlines()[16:]
    cnt = 1
    for line in data_lines:
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].field = line[0]
        sources[-1].freq = float(line[1]) # MHz
        sources[-1].name = 'VWT2014_'+line[2]+'_'+line[0]+'_'+line[1]
        sources[-1].ra = float(line[3]) # deg
        sources[-1].e_ra = float(line[4])/3600.0 # deg (1sigma)
        sources[-1].dec = float(line[5]) # deg
        sources[-1].e_dec = float(line[6])/3600.0 # deg (1sigma)
        sources[-1].total_flux = float(line[7]) # mJy
        sources[-1].e_total_flux = float(line[8]) # mJy (1sigma)
        # These are not given in the catalog
        #   - setting to approx half restoring beam size (table 2 in paper):
        sources[-1].maj = 10. # arcsec
        sources[-1].min = 10. # arcsec
        sources[-1].pa = 0. # deg
        sources[-1].catalog_id = catalog_dict['lofar_lba']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} LOFAR_LBA/VWT2014 sources to '
                      'lofar_lba_psql.txt'.format(len(sources)))
    return sources


def read_lofar_hba(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Nord et al. Galactic Center survey catalog and writes
    them into a file in the same directory called lofar_hba_psql.txt
    if the file does not already exist.

    Telescope/frequency: LOFAR 150 MHz

    Spatial resolution: 6''

    """
    catalog_dict['lofar_hba']['telescope'] = 'LOFAR'
    catalog_dict['lofar_hba']['frequency'] = 150
    catalog_dict['lofar_hba']['resolution'] = 6.
    catalog_dict['lofar_hba']['reference'] = 'Williams et al. (2016)'
    
    psqlf = os.path.join(catalogdir, 'lofar_hba_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir, 'LOFAR_HBA.dat')
    sources = []
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        bflag1 = int(line[13])
        eflag = int(line[14])
        bflag2 = int(line[15])
        aflag = int(line[16])
        if bflag1 == 0 and bflag2 == 0 and eflag == 0 and aflag == 0:
            sources.append(CatalogSource())
            sources[-1].id = cnt
            cnt += 1
            sources[-1].name = 'WVR2016_'+line[0]
            sources[-1].ra = float(line[1]) # deg
            sources[-1].e_ra = float(line[2]) # deg 
            if sources[-1].e_ra < 1e-5: sources[-1].e_ra = 1e-5
            sources[-1].dec = float(line[3]) # deg
            sources[-1].e_dec = float(line[4]) # deg 
            if sources[-1].e_dec < 1e-5: sources[-1].e_dec = 1e-5
            sources[-1].otal_flux = float(line[5]) # mJy
            sources[-1].e_total_flux = float(line[6]) # mJy 
            sources[-1].peak_flux = float(line[7]) # mJy/bm
            sources[-1].e_peak_flux = float(line[8]) # mJy/bm
            sources[-1].rms = float(line[10]) # mJy/bm
            # These are not given in the catalog
            #  - setting to approx half restoring beam size:
            sources[-1].maj = 3.5 # arcsec
            sources[-1].min = 3.5 # arcsec
            sources[-1].pa = 0. # deg
            sources[-1].catalog_id = catalog_dict['lofar_hba']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} LOFAR_HBA/WVR2016 sources to '
                      'lofar_hba_psql.txt'.format(len(sources)))
    return sources


def read_nrl_nvss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Nord et al. Galactic Center survey catalog and writes
    them into a file in the same directory called nrl_nvss_psql.txt
    if the file does not already exist.

    NOTE: This is our version of the NVSS catalog after running
    PyBDSF on all NVSS images.

    Telescope/frequency: VLA 1.4 GHz

    Spatial resolution: 45''

    """
    catalog_dict['nrl_nvss']['telescope'] = 'VLA'
    catalog_dict['nrl_nvss']['frequency'] = 1400
    catalog_dict['nrl_nvss']['resolution'] = 45.
    catalog_dict['nrl_nvss']['reference'] = ''
    
    psqlf = os.path.join(catalogdir, 'nrl_nvss_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir,
                         'NRL-NVSS_uniquecatalog_filtered_5.68.FINAL')
    sources = []
    fread = open(fname, 'r')
    # Skip header
    line = fread.readline()
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        sources.append(CatalogSource())
        line = line.split()
        sources[-1].id = cnt
        sources[-1].name = 'NRLNVSS_%s' % line[0]
        cnt+=1
        sources[-1].ra = float(line[1])  # deg
        sources[-1].dec = float(line[2]) # deg
        sources[-1].e_ra = float(line[3]) # deg
        sources[-1].e_dec = float(line[4]) # deg
        sources[-1].total_flux = float(line[5]) # mJy
        sources[-1].e_total_flux = float(line[6]) # mJy
        sources[-1].peak_flux = float(line[7]) # mJy/bm
        sources[-1].e_peak_flux = float(line[8]) # mJy/bm
        sources[-1].maj = float(line[9])*3600. # arcsec
        sources[-1].min = float(line[10])*3600. # arcsec
        sources[-1].pa = float(line[11]) # deg
        sources[-1].field = line[0][:8]
        sources[-1].catalog_id = catalog_dict['nrl_nvss']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} NRL-NVSS sources to nrl_nvss_psql.txt'.
                      format(len(sources)))
    return sources


def read_atlas(return_sources=False):
    """Generates a list of CatalogSource objects from
    the ATLAS survey catalog and writes them into a file in the 
    same directory called atlas_psql.txt if the file does
    not already exist.

    Telescope/frequency: ATLAS 1.4 GHz

    Spatial resolution: 11'' 
    
    """
    catalog_dict['atlas']['telescope'] = 'ATLAS'
    catalog_dict['atlas']['frequency'] = 1400
    catalog_dict['atlas']['resolution'] = 11.
    catalog_dict['atlas']['reference'] = 'Franzen et al. (2015)'
    
    psqlf = os.path.join(catalogdir, 'atlas_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    
    sources = []
    fname = os.path.join(catalogdir, 'atlasdr3.dat.txt')
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        if line[0] == 'E' or line[0] == 'C':
            sources.append(CatalogSource())
            line = line.split('|')
            sources[-1].id = cnt
            sources[-1].name = line[0][8:]
            cnt += 1
            h = int(line[1].split()[0])
            m = int(line[1].split()[1])
            s = float(line[1].split()[2])
            sources[-1].ra = 15. * dms2deg(h, m, s) # deg
            d = int(line[1].split()[3])
            m = int(line[1].split()[4])
            s = float(line[1].split()[5])
            sources[-1].dec = dms2deg(d, m, s) # deg
            sources[-1].total_flux = float(line[14]) # mJy
            sources[-1].e_total_flux = float(line[15]) # mJy
            sources[-1].peak_flux = float(line[12]) # mJy/beam
            sources[-1].e_peak_flux = float(line[13]) # mJy/beam
            sources[-1].rms = float(line[5]) # mJy/beam
            sources[-1].e_ra = float(line[3])/3600.0 # deg
            sources[-1].e_dec = float(line[4])/3600.0 # deg
            #sources[-1].maj = float(line[9]) # arcsec
            #sources[-1].min = float(line[10]) # arcsec
            #sources[-1].pa = 0.0
            sources[-1].field = line[21][:-1]
            sources[-1].catalog_id = catalog_dict['atlas']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} ATLAS sources to atlas_psql.txt'.format(
        len(sources)))
    return sources


def read_lotssdr1(return_sources=False):
    """Generates a list of CatalogSource objects from
    the LOFAR 2m Sky Survey: LoTSS DR1, catalog and writes
    them into a file in the same directory called lotssdr1_psql.txt
    if the file does not already exist.

    Telescope/frequency: LOFAR 150 MHz

    Spatial resolution: 6''

    """
    catalog_dict['lotssdr1']['telescope'] = 'LOFAR'
    catalog_dict['lotssdr1']['frequency'] = 150
    catalog_dict['lotssdr1']['resolution'] = 6.
    catalog_dict['lotssdr1']['reference'] = 'Shimwell et al. (2019)'
    
    psqlf = os.path.join(catalogdir, 'lotssdr1_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir, 'lotssdr1.dat')
    sources = []
    fread = open(fname, 'r')
    cnt = 1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        sources[-1].ra = float(line[1]) # deg
        sources[-1].e_ra = float(line[2])# deg (total 1sigma error)
        sources[-1].dec = float(line[3]) # deg
        sources[-1].e_dec = float(line[4])# deg (total 1 sigma error)
        sources[-1].total_flux = float(line[5]) # mJy
        sources[-1].e_total_flux = float(line[6]) # mJy
        sources[-1].peak_flux = float(line[7]) # mJy/bm
        sources[-1].e_peak_flux = float(line[8]) # mJy/bm
        sources[-1].maj   = float(line[9]) # arcsec
        sources[-1].e_maj = float(line[10]) # arcsec
        sources[-1].min   = float(line[11]) # arcsec
        sources[-1].e_min = float(line[12]) # arcsec
        sources[-1].pa    = float(line[13]) # deg
        sources[-1].e_pa  = float(line[14]) # deg        
        sources[-1].rms = float(line[15]) # mJy/bm
        sources[-1].field = line[16]
        sources[-1].catalog_id = catalog_dict['lotssdr1']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} LoTSS DR1 sources to lotssdr1_psql.txt'.format(
        len(sources)))
    return sources


def read_racs(return_sources=False):
    """Generates a list of CatalogSource objects from
    the Rapid ASKAP Continuum Survey: RACS, catalog and writes
    them into a file in the same directory called racs_psql.txt
    if the file does not already exist.

    Telescope/frequency: ASKAP 887.5 MHz

    Spatial resolution: 25''

    """
    catalog_dict['racs']['telescope'] = 'ASKAP'
    catalog_dict['racs']['frequency'] = 887.5
    catalog_dict['racs']['resolution'] = 25.0
    catalog_dict['racs']['reference'] = 'Hale et al. (2021)'
    
    psqlf = os.path.join(catalogdir, 'racs_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    #fname = os.path.join(catalogdir, 'lotssdr1.dat')
    fname = '/data3/vpipe/vlite-emil/data/RACS/AS110_Derived_Catalogue_racs_dr1_sources_galacticcut_v2021_08_v01_4551.csv'
    sources = []
    fread = open(fname, 'r')
    #read header
    line = fread.readline()
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split(',')
        sources.append(CatalogSource())
        sources[-1].id = int(line[0])
        tmp = line[2].split()
        sources[-1].name = tmp[0]+'_'+tmp[1]
        sources[-1].ra = float(line[8]) # deg
        sources[-1].e_ra = float(line[10])/3600 # deg 
        sources[-1].dec = float(line[9]) # deg
        sources[-1].e_dec = float(line[11])/3600 # deg 
        sources[-1].total_flux = float(line[12]) # mJy
        sources[-1].e_total_flux = float(line[14]) # mJy
        sources[-1].peak_flux = float(line[15]) # mJy/bm
        sources[-1].e_peak_flux = float(line[16]) # mJy/bm
        sources[-1].maj   = float(line[17]) # arcsec
        sources[-1].e_maj = float(line[18]) # arcsec
        sources[-1].min   = float(line[19]) # arcsec
        sources[-1].e_min = float(line[20]) # arcsec
        sources[-1].pa    = float(line[21]) # deg
        sources[-1].e_pa  = float(line[22]) # deg        
        sources[-1].rms = float(line[31]) # mJy/bm
        sources[-1].field = line[4]
        sources[-1].catalog_id = catalog_dict['racs']['id']
        sources[-1].pt_like = False
        frat = sources[-1].total_flux / sources[-1].peak_flux
        snr  = sources[-1].peak_flux / sources[-1].rms
        if frat < (1.025 +0.69*pow(snr,-0.62)):
            sources[-1].pt_like = True
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} RACS sources to racs_psql.txt'.format(
        len(sources)))
    return sources


def read_psr(return_sources=False):
    """Generates a list of CatalogSource objects from
    the ATNF pulsar catalog and writes them into a file in the same 
    directory called psr_psql.txt if the file does not already exist.

    Telescope/frequency: Many, set to: ATNF. Will include 1400 MHz fluxes,
      if given

    Spatial resolution: No fixed resolution, set to: 0''

    """
    catalog_dict['psr']['telescope'] = 'ATNF'
    catalog_dict['psr']['frequency'] = 1400.
    catalog_dict['psr']['resolution'] = 0.
    catalog_dict['psr']['reference'] = 'version 1.65'
    
    psqlf = os.path.join(catalogdir, 'psr_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass

    fname = os.path.join(catalogdir, 'psrcat_short.txt')
    sources = []
    fread = open(fname, 'r')
    #read headers
    line = fread.readline()
    line = fread.readline()
    line = fread.readline()
    line = fread.readline()
    cnt=1
    while 1:
        line = fread.readline()
        if not line: break
        line = line.split()
        sources.append(CatalogSource())
        sources[-1].id = cnt
        cnt += 1
        sources[-1].name = line[0]
        sources[-1].ra = float(line[1]) # deg
        sources[-1].e_ra = float(line[2]) # deg 
        sources[-1].dec = float(line[3]) # deg
        sources[-1].e_dec = float(line[4]) # deg
        try:
            sources[-1].total_flux = float(line[5]) # mJy
            try:
                sources[-1].e_total_flux = float(line[6]) # mJy
            except:
                sources[-1].e_total_flux = None
        except:
             sources[-1].total_flux = None
             sources[-1].e_total_flux = None
        sources[-1].catalog_id = catalog_dict['psr']['id']
    fread.close()
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i %s None None None\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id, src.pt_like))
    catio_logger.info(' -- wrote {} ATNF pulsars sources to psr_psql.txt'.format(
        len(sources)))
    return sources

def read_vlass1(return_sources=False):
    """Generates a list of CatalogSource objects from
    the VLASS epoch 1 QuickLook catalog and writes them into a file in the 
    same directory called vlass1_psql.txt if the file does
    not already exist.

    Telescope/frequency: VLA 3000 MHz

    Spatial resolution: 2.5" 

    """
    catalog_dict['vlass1']['telescope'] = 'VLA'
    catalog_dict['vlass1']['frequency'] = 3000.
    catalog_dict['vlass1']['resolution'] = 2.5
    catalog_dict['vlass1']['reference'] = 'cirada.ca/vlasscatalogueql0'
    
    psqlf = os.path.join(catalogdir, 'vlass1_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        #read file
        sources = pd.read_csv('/roadraid2/emilp/VLASS/epoch1/CIRADA_VLASS1QL_table1_components_v2.csv.gz')
        #output file
        psqlf = 'vlass1_psql.txt'
        #loop over srcs & write to file
        with open(psqlf, 'w') as fwrite:
            for i,row in sources.iterrows():
                fwrite.write('%d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %i %s %s %s %s\n' % (i+1,row['Component_name'],row['RA'],row['E_RA'],row['DEC'],row['E_DEC'],row['Total_flux'],row['E_Total_flux'],row['Peak_flux'],row['E_Peak_flux'],row['Maj'],row['E_Maj'],row['Min'],row['E_Min'],row['PA'],row['E_PA'],row['Isl_rms'],row['Tile'],25,'None',row['Duplicate_flag'],row['Quality_flag'],row['P_sidelobe']))
        catio_logger.info(' -- wrote {} sources to vlass1_psql.txt'.format(len(sources)))
        return sources
    
        
def read_vlass2(return_sources=False):
    """Generates a list of CatalogSource objects from
    the VLASS epoch 2 QuickLook catalog and writes them into a file in the 
    same directory called vlass1_psql.txt if the file does
    not already exist.

    Telescope/frequency: VLA 3000 MHz

    Spatial resolution: 2.5" 

    """
    catalog_dict['vlass2']['telescope'] = 'VLA'
    catalog_dict['vlass2']['frequency'] = 3000.
    catalog_dict['vlass2']['resolution'] = 2.5
    catalog_dict['vlass2']['reference'] = 'cirada.ca/vlasscatalogueql0'
    
    psqlf = os.path.join(catalogdir, 'vlass2_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        #read file
        sources=pd.read_csv('/roadraid2/emilp/VLASS/epoch2/CIRADA_VLASS2QLv1_table1_components.csv.gz')
        #output file
        psqlf = 'vlass2_psql.txt'
        #loop over srcs & write to file
        with open(psqlf, 'w') as fwrite:
            for i,row in sources.iterrows():
                fwrite.write('%d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %i %s %s %s %s\n' % (i+1,row['Component_name'],row['RA'],row['E_RA'],row['DEC'],row['E_DEC'],row['Total_flux'],row['E_Total_flux'],row['Peak_flux'],row['E_Peak_flux'],row['Maj'],row['E_Maj'],row['Min'],row['E_Min'],row['PA'],row['E_PA'],row['Isl_rms'],row['Tile'],25,'None',row['Duplicate_flag'],row['Quality_flag'],'None'))
        catio_logger.info(' -- wrote {} sources to vlass2_psql.txt'.format(len(sources)))
        return sources


'''
def read_vlass3(return_sources=False):
    """Generates a list of CatalogSource objects from
    the VLASS epoch 2 QuickLook catalog and writes them into a file in the 
    same directory called vlass1_psql.txt if the file does
    not already exist.

    Telescope/frequency: VLA 3000 MHz

    Spatial resolution: 2.5" 

    """
    catalog_dict['vlass3']['telescope'] = 'VLA'
    catalog_dict['vlass3']['frequency'] = 3000.
    catalog_dict['vlass3']['resolution'] = 2.5
    catalog_dict['vlass3']['reference'] = 'cirada.ca/vlasscatalogueql0'
    
    psqlf = os.path.join(catalogdir, 'vlass3_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        #read file
        sources=pd.read_csv('/roadraid2/emilp/VLASS/epoch3/CIRADA_VLASS3QLv1_table1_components.csv.gz')
        #output file
        psqlf = 'vlass3_psql.txt'
        #loop over srcs & write to file
        with open(psqlf, 'w') as fwrite:
            for i,row in sources.iterrows():
                fwrite.write('%d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %i %s %s %s %s\n' % (i+1,row['Component_name'],row['RA'],row['E_RA'],row['DEC'],row['E_DEC'],row['Total_flux'],row['E_Total_flux'],row['Peak_flux'],row['E_Peak_flux'],row['Maj'],row['E_Maj'],row['Min'],row['E_Min'],row['PA'],row['E_PA'],row['Isl_rms'],row['Tile'],25,'None',row['Duplicate_flag'],row['Quality_flag'],'None'))
        catio_logger.info(' -- wrote {} sources to vlass3_psql.txt'.format(len(sources)))
        return sources
'''
