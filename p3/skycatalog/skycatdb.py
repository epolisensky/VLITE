"""Tables corresponding to each sky survey catalog are 
stored in a separate schema called 'skycat' within the 
same `PostgreSQL` database as the VLITE data. This
module contains the functions to create the schema
and tables.

Currently included catalogs:
1. FIRST - VLA 1.4 GHz, 5" res.
2. GLEAM - MWA 74-231 MHz, 100" res.
3. NVSS - VLA 1.4 GHz, 45" res.
4. SUMSS - MOST 843 MHz, 45" res.
5. TGSS - GMRT 150 MHz, 25" res.
6. WENSS - WSRT 325 MHz, 54" res.

"""
import os
import psycopg2
from psycopg2 import sql
import catalogio


def index(tblname, conn):
    """Indexes and clusters the table using Q3C.

    Parameters
    ----------
    tblname : str
        Name of the database table.
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    """
    cur = conn.cursor()

    # Just in case...
    tblname = tblname.lower()
        
    cur.execute(psycopg2.sql.SQL('''CREATE INDEX ON skycat.{} 
        (q3c_ang2ipix(ra, dec))''').format(psycopg2.sql.Identifier(tblname)))
    cur.execute(psycopg2.sql.SQL('CLUSTER {} ON skycat.{}').format(
        psycopg2.sql.Identifier(tblname + '_q3c_ang2ipix_idx'),
        psycopg2.sql.Identifier(tblname)))
    cur.execute(psycopg2.sql.SQL('ANALYZE skycat.{}').format(
        psycopg2.sql.Identifier(tblname)))

    conn.commit()
    cur.close()


def add_table(tblname, conn):
    """Copies sky survey sources into a new table
    in the 'skycat' schema from a text file. If the
    text file does not yet exist, the sources are
    read from the original catalog, initialized as
    CatalogSource objects, and written to the text file.

    Parameters
    ----------
    tblname : str
        Name of the table/catalog from which the source
        originated.
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.
    """
    cur = conn.cursor()

    # force to all lowercase, just in case
    tblname = tblname.lower()
   
    # Create the table
    cur.execute(psycopg2.sql.SQL(
        '''
        CREATE TABLE IF NOT EXISTS skycat.{} (
            id INTEGER NOT NULL,
            name TEXT,
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
            catalog_id INTEGER,
            PRIMARY KEY (id),
            FOREIGN KEY (catalog_id)
              REFERENCES skycat.catalogs (id)
              ON UPDATE CASCADE
              ON DELETE CASCADE
        );
        ''').format(psycopg2.sql.Identifier(tblname)))

    # Read the sky catalog
    print('Adding sources from {}'.format(tblname))
    if tblname == 'tgss':
        sources = catalogio.read_tgss()
    elif tblname == 'nvss':
        sources = catalogio.read_nvss()
    elif tblname == 'first':
        sources = catalogio.read_first()
    elif tblname == 'sumss':
        sources = catalogio.read_sumss()
    elif tblname == 'wenss':
        sources = catalogio.read_wenss()
    elif tblname == 'gleam':
        sources = catalogio.read_gleam()
    elif tblname == 'cosmos':
        sources = catalogio.read_cosmos()
    elif tblname == 'vlssr':
        sources = catalogio.read_vlssr()
    elif tblname == 'txs':
        sources = catalogio.read_txs()
    elif tblname == 'sevenc':
        sources = catalogio.read_sevenc()
    elif tblname == 'gpsr5':
        sources = catalogio.read_gpsr5()
    elif tblname == 'gpsr1':
        sources = catalogio.read_gpsr1()
    elif tblname == 'nordgc':
        sources = catalogio.read_nordgc()
    else:
        print('\nERROR: No function to read catalog {}.'.format(tblname))

    # Add catalog to catalogs table
    sql = '''INSERT INTO skycat.catalogs (id, name) VALUES (%s, %s)
        ON CONFLICT (id) DO NOTHING'''
    cur.execute(sql, (catalogio.catid_dict[tblname], tblname))

    # This way took too long (~11 min)
    """
    # Add sources to table
    for src in sources:
        src.catalog_id = catid
        sql = '''INSERT INTO skycat.{} (
            name, ra, e_ra, dec, e_dec, total_flux, e_total_flux,
            peak_flux, e_peak_flux, maj, e_maj, min, e_min, pa, e_pa,
            rms, field, catalog_id) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s)'''
        vals = (src.name, src.ra, src.e_ra, src.dec, src.e_dec,
                src.total_flux, src.e_total_flux, src.peak_flux,
                src.e_peak_flux, src.maj, src.e_maj, src.min, src.e_min,
                src.pa, src.e_pa, src.rms, src.field, src.catalog_id)
        cur.execute(psycopg2.sql.SQL(sql).format(
            psycopg2.sql.Identifier(tblname)), vals)
    """

    # This is waaaaaaay faster (~10x)
    fname = tblname + '_psql.txt'
    psqlf = os.path.join(catalogio.catalogdir, fname)
    fullname = 'skycat.' + tblname
    with open(psqlf, 'r') as f:
        cur.copy_from(f, fullname, sep=' ', null='None')

    conn.commit()
    cur.close()
    

def create(conn):
    """Creates the 'skycat' schema and catalogs 
    table within that schema. Calls the functions
    `add_table` and `index` for each catalog in the
    list from `catalogio.catalog_list`.

    Parameters
    ----------
    conn : psycopg2.extensions.connect instance
        The `PostgreSQL` database connection object.    
    """
    cur = conn.cursor()
    sql = (
        '''
        CREATE SCHEMA IF NOT EXISTS skycat;
        
        CREATE TABLE IF NOT EXISTS skycat.catalogs (
            id INTEGER NOT NULL,
            name TEXT,
            PRIMARY KEY (id)
        );
        ''')
    cur.execute(sql)
    cur.close

    tables = catalogio.catalog_list
    for table in tables:
        add_table(table, conn)
        index(table, conn)
