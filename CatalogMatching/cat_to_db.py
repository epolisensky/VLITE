"""cat_to_db.py reads in the sky survey catalogs 
using the catalogio.py module and stores them
in separate but identical tables in a database
called SkyCatalogs.sqlite.

Currently included catalogs:
1. TGSS - GMRT 150 MHz, 25" res.
2. NVSS - VLA 1.4 GHz, 45" res.
3. FIRST - VLA 1.4 GHz, 5" res.
4. SUMSS - MOST 843 MHz, 45" res.
5. WENSS - WSRT 325 MHz, 54" res.
6. GLEAM - MWA 74-231 MHz, 100" res."""


import sqlite3
import sys
import os
import pandas as pd
import catalogio


def fill_null(catalog, attrs):
    for attr in attrs:
        s = pd.Series([getattr(src, attr) for src in catalog])
        s.loc[s.notnull() == False] = s.median()
        for idx in range(len(catalog)):
            setattr(catalog[idx], attr, s[idx])
    return catalog


# Read sky catalogs
print 'Reading TGSS...'
tgss = catalogio.readTGSS()

print 'Reading NVSS...'
nvss = fill_null(catalogio.readNVSS(), ['e_ra', 'e_dec'])

print 'Reading FIRST...'
first = catalogio.readFIRST()

print 'Reading SUMSS...'
sumss = catalogio.readSUMSS()

print 'Reading WENSS...'
wenss = catalogio.readWENSS()

print 'Reading GLEAM...'
gleam = fill_null(catalogio.readGLEAM(), ['e_ra', 'e_dec'])

tables = ['TGSS', 'NVSS', 'FIRST', 'SUMSS', 'WENSS', 'GLEAM']
catalogs = [tgss, nvss, first, sumss, wenss, gleam]

catid_dict = {}
idnum = 1
for table in sorted(tables):
    catid_dict[table] = idnum
    idnum += 1

dbname = os.path.join(catalogio.catalogdir, 'SkyCatalogs.sqlite')
conn = sqlite3.connect(dbname)
cur = conn.cursor()

cur.executescript('''
DROP TABLE IF EXISTS CatalogID;

CREATE TABLE CatalogID (
    id INTEGER NOT NULL UNIQUE PRIMARY KEY,
    name TEXT UNIQUE
    );
''')
for catalog in catid_dict.keys():
    cur.execute('INSERT INTO CatalogID (id, name) VALUES (?, ?)',
                (catid_dict[catalog], catalog))

for table in tables:
    cur.execute('DROP TABLE IF EXISTS %s' % table)

    cur.execute('''CREATE TABLE %s (
        id INTEGER NOT NULL UNIQUE PRIMARY KEY,
        name TEXT UNIQUE,
        ra REAL,
        e_ra REAL,
        dec REAL,
        e_dec REAL,
        total_flux REAL,
        e_total_flux REAL,
        peak_flux REAL,
        e_peak_flux REAL,
        maj REAL,
        e_maj REAL,
        min REAL,
        e_min REAL,
        pa REAL,
        e_pa REAL,
        rms REAL,
        field TEXT,
        catalog_id INTEGER
    )''' % table)


for src in nvss:
    cur.execute('''INSERT OR IGNORE INTO NVSS (
        name, ra, e_ra, dec, e_dec, total_flux, e_total_flux,
        peak_flux, e_peak_flux, maj, e_maj, min, e_min, pa, e_pa,
        rms, field, catalog_id) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (src.name, src.ra, src.e_ra, src.dec, src.e_dec,
                 src.total_flux, src.e_total_flux, src.peak_flux,
                 src.e_peak_flux, src.maj, src.e_maj, src.min,
                 src.e_min, src.pa, src.e_pa, src.rms, src.field,
                 catid_dict['NVSS']))

for src in tgss:
    cur.execute('''INSERT OR IGNORE INTO TGSS (
        name, ra, e_ra, dec, e_dec, total_flux, e_total_flux,
        peak_flux, e_peak_flux, maj, e_maj, min, e_min, pa, e_pa,
        rms, field, catalog_id) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (src.name, src.ra, src.e_ra, src.dec, src.e_dec,
                 src.total_flux, src.e_total_flux, src.peak_flux,
                 src.e_peak_flux, src.maj, src.e_maj, src.min,
                 src.e_min, src.pa, src.e_pa, src.rms, src.field,
                 catid_dict['TGSS']))

for src in first:
    cur.execute('''INSERT OR IGNORE INTO FIRST (
        name, ra, e_ra, dec, e_dec, total_flux, e_total_flux,
        peak_flux, e_peak_flux, maj, e_maj, min, e_min, pa, e_pa,
        rms, field, catalog_id) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (src.name, src.ra, src.e_ra, src.dec, src.e_dec,
                 src.total_flux, src.e_total_flux, src.peak_flux,
                 src.e_peak_flux, src.maj, src.e_maj, src.min,
                 src.e_min, src.pa, src.e_pa, src.rms, src.field,
                 catid_dict['FIRST']))

for src in sumss:
    cur.execute('''INSERT OR IGNORE INTO SUMSS (
        name, ra, e_ra, dec, e_dec, total_flux, e_total_flux,
        peak_flux, e_peak_flux, maj, e_maj, min, e_min, pa, e_pa,
        rms, field, catalog_id) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (src.name, src.ra, src.e_ra, src.dec, src.e_dec,
                 src.total_flux, src.e_total_flux, src.peak_flux,
                 src.e_peak_flux, src.maj, src.e_maj, src.min,
                 src.e_min, src.pa, src.e_pa, src.rms, src.field,
                 catid_dict['SUMSS']))

for src in wenss:
    cur.execute('''INSERT OR IGNORE INTO WENSS (
        name, ra, e_ra, dec, e_dec, total_flux, e_total_flux,
        peak_flux, e_peak_flux, maj, e_maj, min, e_min, pa, e_pa,
        rms, field, catalog_id) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (src.name, src.ra, src.e_ra, src.dec, src.e_dec,
                 src.total_flux, src.e_total_flux, src.peak_flux,
                 src.e_peak_flux, src.maj, src.e_maj, src.min,
                 src.e_min, src.pa, src.e_pa, src.rms, src.field,
                 catid_dict['WENSS']))

for src in gleam:
    cur.execute('''INSERT OR IGNORE INTO GLEAM (
        name, ra, e_ra, dec, e_dec, total_flux, e_total_flux,
        peak_flux, e_peak_flux, maj, e_maj, min, e_min, pa, e_pa,
        rms, field, catalog_id) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (src.name, src.ra, src.e_ra, src.dec, src.e_dec,
                 src.total_flux, src.e_total_flux, src.peak_flux,
                 src.e_peak_flux, src.maj, src.e_maj, src.min,
                 src.e_min, src.pa, src.e_pa, src.rms, src.field,
                 catid_dict['GLEAM']))

conn.commit()
cur.close()
