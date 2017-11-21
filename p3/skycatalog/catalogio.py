"""This module contains all the messy I/O functionality
for reading in the different sky survey catalogs.
Each catalog source is initialized as a CatalogSource
object to create a uniform set of properties for each
sky survey's sources. When reading the catalogs for the
first time, sources with the uniform set of attributes
are written to text files called '[catalog_name]_psql.txt'
which is used to insert the sources into the `PostgreSQL`
database in an efficient way. 

"""
import os
import pandas as pd


catalogdir = '/home/erichards/work/data/SkyCatalogs'

# If adding a new catalog, don't forget to add it to this list!
catalog_list = sorted(['tgss', 'nvss', 'first', 'sumss', 'wenss', 'gleam'])

catid_dict = {}
id = 1
for catalog in catalog_list:
    catid_dict[catalog] = id
    id += 1


def dms2deg(d, m, s):
    """Translates coordinates from deg:min:sec to decimal degrees.
    If computing RA, the result needs to be multiplied by 15.

    """
    dec = abs(float(d)) + abs(float(m))/60. + abs(float(s))/3600.
    d = str(d)
    if d.__contains__("-"):
        dec = -dec
    return dec


def set_error(catalog, attrs):
    """Sets missing values for specified attributes to
    the median value of the distribution.

    """
    for attr in attrs:
        s = pd.Series([getattr(src, attr) for src in catalog])
        s.loc[s.notnull() == False] = s.median()
        for idx in range(len(catalog)):
            setattr(catalog[idx], attr, s[idx])
    return catalog


class CatalogSource(object):
    """Class for the sky survey sources. Their
    catalog origin is identified by the catalog_id
    attribute.
    
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


def read_tgss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the TGSS survey catalog and writes them into a file in the 
    same directory called tgss_psql.txt if the file does
    not already exist.

    Telescope/frequency: GMRT 150 MHz
    Spatial resolution: 25'' 

    """
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
        sources[-1].catalog_id = catid_dict['tgss']
    print 'Read %d TGSS sources' % len(sources)
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%i %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id))
    return sources


def read_first(return_sources=False):
    """Generates a list of CatalogSource objects from
    the FIRST survey catalog and writes them into a file in the 
    same directory called first_psql.txt if the file does
    not already exist.

    Telescope/frequency: VLA 1.4 GHz
    Spatial resolution: 5'' 

    """
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
                sources[-1].ra = float(15*dms2deg(line[0],
                                                  line[1], line[2])) # deg
                sources[-1].dec = float(dms2deg(line[3],
                                                line[4], line[5])) # deg
                sources[-1].prob = float(line[6])
                sources[-1].peak_flux = float(line[7]) # mJy/beam
                sources[-1].total_flux = float(line[8]) # mJy
                sources[-1].rms = float(line[9]) # mJy/beam
                sources[-1].maj = float(line[10]) # arcsec
                sources[-1].min = float(line[11]) # arcsec
                sources[-1].pa = float(line[12])
                # HACK! setting FIRST positional uncertainty to 1 arcsec
                # actual uncertainties need to be calculated, see FIRST
                # website or paper
                # sources[-1].e_ra  = 1.0/3600.0 # deg
                # sources[-1].e_dec = 1.0/3600.0 # deg
                # HACK: FIRST positional uncertainty calculations given
                # for FITTED bmaj & bmin. 
                #  I'll assign these to e_ra & e_dec for lack of better method
                snr = (sources[-1].peak_flux-0.25) / sources[-1].rms
                dfMaj = float(line[13]) * ((1.0/snr) + 0.05) / 3600.0 # deg
                dfMin = float(line[14]) * ((1.0/snr) + 0.05) / 3600.0 # deg
                # fpa = float(line[15]) # deg
                sources[-1].e_ra = dfMin # deg
                sources[-1].e_dec = dfMaj # deg
                sources[-1].catalog_id = catid_dict['first']
    fread.close()
    print 'Read %d FIRST sources' % len(sources)
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%i %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id))
    return sources


def read_sumss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the SUMSS survey catalog and writes them into a file in the 
    same directory called sumss_psql.txt if the file does
    not already exist.

    Telescope/frequency: MOST 843 MHz
    Spatial resolution: 45'' 

    """
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
            sources[-1].ra = float(15.0*dms2deg(line[0],
                                                line[1], line[2])) # deg
            sources[-1].dec = float(dms2deg(line[3],
                                            line[4], line[5])) # deg
            sources[-1].e_ra = float(line[6])/3600.0 # deg
            sources[-1].e_dec = float(line[7])/3600.0 # deg
            sources[-1].peak_flux = float(line[8]) # mJy/bm
            sources[-1].e_peak_flux = float(line[9])
            sources[-1].total_flux = float(line[10]) # mJy
            sources[-1].e_total_flux = float(line[11])
            sources[-1].maj = float(line[12]) # arcsec
            sources[-1].min = float(line[13]) # arcsec
            sources[-1].pa = float(line[14])
            sources[-1].catalog_id = catid_dict['sumss']
    fread.close()
    print 'Read %d SUMSS sources' % len(sources)
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%i %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id))
    return sources


def read_wenss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the WENSS survey catalog and writes them into a file in the 
    same directory called wenss_psql.txt if the file does
    not already exist.

    Telescope/frequency: WSRT 325 MHz
    Spatial resolution: 54'' 
    
    """
    psqlf = os.path.join(catalogdir, 'wenss_psql.txt')
    if os.path.isfile(psqlf):
        if not return_sources:
            return
        else:
            pass
    else:
        pass
    sources = []
    fname = os.path.join(catalogdir, 'WENS.COMPLETE.txt')
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
            sources[-1].catalog_id = catid_dict['wenss']
    fread.close()
    print 'Read %d WENSS sources' % len(sources)
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%i %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id))
    return sources


def read_nvss(return_sources=False):
    """Generates a list of CatalogSource objects from
    the NVSS survey catalog and writes them into a file in the 
    same directory called nvss_psql.txt if the file does
    not already exist.

    Telescope/frequency: VLA 1.4 GHz
    Spatial resolution: 45'' 

    """
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
                sources[-1].ra = float(15*dms2deg(line[0:2], line[3:5],
                                                  line[7:11]))
            else:
                sources[-1].ra = float(15*dms2deg(line[0:2], line[3:5],
                                                  line[6:11]))
            if line[19] == ' ': # dec
                sources[-1].dec = float(dms2deg(line[12:15], line[16:18],
                                                line[20:23]))
            else:
                sources[-1].dec = float(dms2deg(line[12:15], line[16:18],
                                                line[19:23]))
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
            sources[-1].catalog_id = catid_dict['nvss']
    fread.close()
    print 'Read %d NVSS sources' % len(sources)
    sources = set_error(sources, ['e_ra', 'e_dec'])
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%s %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %s\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id))
    return sources


def read_gleam(return_sources=False):
    """Generates a list of CatalogSource objects from
    the GLEAM survey catalog and writes them into a file in the 
    same directory called gleam_psql.txt if the file does
    not already exist.

    Telescope/frequency: MWA 74-231 MHz
    Spatial resolution: ~100'' 

    """
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
        sources[-1].catalog_id = catid_dict['gleam']
    print 'Read %d GLEAM sources' % len(sources)
    sources = set_error(sources, ['e_ra', 'e_dec'])
    with open(psqlf, 'w') as fwrite:
        for src in sources:
            fwrite.write('%i %s %s %s %s %s %s %s %s %s %s %s %s '
                         '%s %s %s %s %s %i\n' % (
                             src.id, src.name, src.ra, src.e_ra, src.dec,
                             src.e_dec, src.total_flux, src.e_total_flux,
                             src.peak_flux, src.e_peak_flux, src.maj,
                             src.e_maj, src.min, src.e_min, src.pa, src.e_pa,
                             src.rms, src.field, src.catalog_id))
    return sources
