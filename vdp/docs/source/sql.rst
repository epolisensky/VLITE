.. _postgresql:

PostgreSQL Primer
=================
While a full description of PostgreSQL is way beyond the
scope of this documentation, included here is a short
primer on PostgreSQL, its interactive terminal, and
example SQL queries. Consult the `PostgreSQL documentation
<https://www.postgresql.org/docs/10/static/index.html>`_
for your remaining SQL needs.

See "postgresql_install_notes.txt" on roadrunner in
the vpipe home directory for installation instructions.

.. _server:

Managing the Database Server
----------------------------
A database storage area on disk called a database cluster
was initialized when PostgreSQL was installed on roadrunner.
A database cluster is just "a collection of databases that is
managed by a single instance of a running database server,"
according to the PostgreSQL docs. This is essentially where
all data stored in the database lives and the location is
referenced by the environment variable PGDATA. The data directory
is typically called "data" and must be owned by the
postgres user.

Starting and Stopping the Server
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The server application ``pg_ctl`` can be used to easily start
and stop the database server. The server should be shutdown
any time its host will be shutdown. It will need to be re-started
after every shutdown. You must be logged in as superuser
postgres to run the pg_ctl commands::

  $ su - postgres

To check if the server is already running::
  
  $ pg_ctl status

To start the server::
  
  $ pg_ctl start -l logfile

This command assumes the PGDATA environment variable has been set
and the location of the pgsql/bin/ directory has been added to the path
for the postgres user. This has already been done on roadrunner.
The database log file will be written in the postgres user home directory
(or whatever directory you happened to be in when you issued the command).
You can specifiy a different directory for the log if you so desire.

To shutdown the server::

  $ pg_ctl stop


.. _psql:

Connecting to a Database
------------------------
For fast and easy access to a database, you can use the
PostgreSQL interactive terminal, ``psql``. ``psql`` allows
you to issue any SQL statement and also comes with its
own set of useful meta-commands. To start an interactive
session, just type "psql" and the name of the database
you wish to connect to while logged in as an existing
database user::

  $ psql vlite16

To quit the interactive session, type "\\q".

.. note:: It is *strongly recommended* that you connect
	  to the database while logged in under your
	  own user account rather than as the vpipe user.
	  The permissions for all users other than vpipe and
	  postgres have been set such that you can query the
	  database but cannot alter it in any way. If you
	  choose to be vpipe, just remember with great
	  power comes great responsibility. Query wisely.

You will need to add the pgsql/bin directory to your user's
path in order to run PostgreSQL's client-side applications
like psql without supplying the full directory path::

  PATH=/usr/local/pgsql/bin:$PATH
  export PATH

Summarized below are some psql meta-commands you might find useful:

\\df[tn]
  List all functions. The t option lists only "trigger"
  functions and the n option lists only "normal"
  functions. With **vdp**, there should be 25 "normal"
  functions created by Q3C and 4 "trigger" functions.

\\di
  List all indexes.

\\dn
  List all the schemas.

\\dt
  List all tables in the "public" schema.

\\du or \\dg
  List all the database users/roles.

\\dx
  List all extensions. You can use this to verify
  that Q3C is enabled.

\\l[+]
  List all databases in the server. If + is appended, then
  the sizes of the databases will also be shown. You
  can also run this command outside the interactive
  session with ``psql -l``.

\\q
  Quit.

There is also the option to connect to a database in a
non-interactive way through Python using the ``psycopg2``
module, which is how **vpd** operates.

Copying Data to a File
----------------------
It is possible to export data from the database by
copying it to a text, csv, or binary file. This can be accomplished
through the frontend (client-side) \copy or backend (server-side)
SQL COPY, the latter being more efficient and preferable
for large amounts of data. The client-side copy accepts both
relative and absolute file paths. Accessibility and privileges
are those of the local user, not the server, and no SQL superuser
privileges are required. The server-side copy accepts only absolute
file paths and the location must be accessible by the postgres user.
Therefore, for large copies, you will need to connect to the database
as the vpipe or postgres user.

To copy an entire table to a text, csv, or binary file:

- client-side:

  .. code-block:: psql

	vlite16=> \copy image to 'images.txt';
	vlite16=> \copy image to '/home/erichards/images.csv' csv;
	vlite16=> \copy image to 'images.bin' binary;

- server-side:

  .. code-block:: sql

	vlite16=# COPY image TO '/home/postgres/images.txt';
	vlite16=# COPY image TO '/home/postgres/images.csv' (FORMAT csv);
	vlite16=# COPY image TO '/home/postgres/images.bin' (FORMAT binary);

To copy the result of a query to a text file:

- client-side:

  .. code-block:: psql

	vlite16=> \copy (SELECT * FROM detected_source WHERE image_id = 1) to 'image1sources.txt';

- server-side:

  .. code-block:: sql

	vlite16=# COPY (SELECT * FROM detected_source WHERE image_id = 1) TO '/home/postgres/image1sources.txt';

.. _queries:

Query Examples
--------------
The following query examples can be used to extract
information from any database created by the VLITE
Database Pipeline while connected through an interactive
``psql`` session.

Basic SQL
^^^^^^^^^
.. highlight:: sql
	       
- Display entire table contents::

    TABLE image;
    
  or::

    SELECT * FROM image;

- Limit to certain number of rows::

    SELECT * FROM image LIMIT 2;

- Order the results by a column::

    SELECT * FROM image ORDER BY id LIMIT 2;

- Only select certain columns::

    SELECT id, filename, config, noise, nsrc FROM image;

- Group results on a single or multiple columns::

    SELECT config, obs_date FROM image GROUP BY (config, obs_date)
    ORDER BY obs_date;
    
    config |  obs_date  
    --------+------------
     C      | 2017-07-25
     C      | 2017-07-29
     C      | 2017-08-01
     C      | 2017-08-08
     B      | 2018-01-07
     A      | 2018-03-26
    (6 rows)

- Select rows where a column is equal to some value::

    SELECT * FROM image WHERE error_id IS NOT NULL;
    SELECT * FROM image WHERE config = 'A';

- Count the number of rows::

    SELECT COUNT(1) FROM image;
    SELECT COUNT(1) FROM detected_source WHERE image_id = 1;

- Select rows where a column is in a range of values::

    SELECT * FROM image WHERE nsrc BETWEEN 400 AND 500;

- Select rows where a column is in a list of values::

    SELECT * FROM detected_source WHERE image_id IN (1, 3, 5);

- Select rows where multiple columns are equal to some value::

    SELECT * FROM detected_source WHERE image_id = 1 AND src_id = 0;
    
  or::

    SELECT * FROM detected_source WHERE (image_id, src_id) = (1, 0);

- Select rows using string pattern matching::

    SELECT * FROM image WHERE filename LIKE '%2017-07%';
    
  The '%' symbols act as wildcards.

- Join tables::

    SELECT ra, dec, maj, min, cf.total_flux, snr
    FROM detected_source ds JOIN corrected_flux cf ON
        (ds.src_id, ds.image_id) = (cf.src_id, cf.image_id)
    WHERE distance_from_center < 1.5;
    
  or::

    SELECT a.ra, a.dec, a.maj, a.min, b.total_flux, b.snr
    FROM detected_source AS a, corrected_flux AS b
    WHERE (a.src_id, a.image_id) = (b.src_id, b.image_id) AND
        b.distance_from_center < 1.5;

- Subqueries (sort of like nested queries)::

    SELECT * FROM assoc_source WHERE ndetect = (
        SELECT MAX(ndetect) FROM assoc_source);
	
    SELECT * FROM image WHERE id IN (
        SELECT image_id FROM vlite_unique WHERE detected);
    

Queries for the VLITE Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Create a list of sources from a single image with corrected fluxes::

    SELECT ds.src_id, ra, e_ra, dec, e_dec, cf.total_flux, cf.e_total_flux,
        cf.peak_flux, cf.e_peak_flux, maj, min, pa, code,
	distance_from_center, snr
    FROM detected_source ds JOIN corrected_flux cf ON
        (ds.src_id, ds.image_id) = (cf.src_id, cf.image_id)
    WHERE ds.image_id = 1;

- Select all sources within 1 degree of a position (cone search)::

    SELECT * FROM assoc_source WHERE q3c_join(
        39.9704166625, -1.576805555, ra, dec, 1.0);
    
  When the number of sources (rows) starts to get in to the millions,
  use this query instead::

    SELECT * FROM assoc_source WHERE q3c_radial_query(
        ra, dec, 39.9704166625, -1.576805555, 1.0);
    
  See the `Q3C page <https://github.com/segasai/q3c>`_ for details
  and more queries.

- Get all individual detections with corrected fluxes for a
  single association::

    SELECT ds.src_id, ds.image_id, ra, e_ra, dec, e_dec, cf.total_flux,
        cf.e_total_flux, cf.peak_flux, cf.e_peak_flux, maj, min, pa,
        distance_from_center, snr
    FROM detected_source ds JOIN corrected_flux cf ON
        (ds.src_id, ds.image_id) = (cf.src_id, cf.image_id)
    WHERE assoc_id = 1482;

  Compare properties of the different images which contain the sources::

    SELECT * FROM image WHERE id IN (
        SELECT image_id FROM detected_source WHERE assoc_id = 1482);

- Get all detections of a source, including missed associations and
  detections at other spatial resolutions by peforming a cone search
  around the source's position::

    SELECT * FROM detected_source WHERE q3c_join(
        (SELECT ra FROM assoc_source WHERE id = 1482),
        (SELECT dec FROM assoc_source WHERE id = 1482),
	ra, dec, (10./3600.));

- Get corrected fluxes and MJD times for all the above detections::

    SELECT mjdtime, cf.total_flux FROM image
    JOIN corrected_flux cf ON
        image.id = cf.image_id
    JOIN detected_source ds ON
        (cf.src_id, cf.image_id) = (ds.src_id, ds.image_id)
    WHERE q3c_join(
        (SELECT ra FROM assoc_source WHERE id = 1482),
        (SELECT dec FROM assoc_source WHERE id = 1482),
	ra, dec, (10./3600.))
    ORDER BY mjdtime;
    
     id  | config |  bmaj   |  bmin   |    mjdtime    |    total_flux    
    -----+--------+---------+---------+---------------+------------------
      20 | C      | 72.3309 | 57.6424 | 57959.8536227 | 1020.23183185697
      22 | C      |  69.074 | 57.8562 | 57959.8762385 | 800.673916800101
      24 | C      | 77.1984 | 58.3155 | 57959.9078009 | 1278.92413095431
      28 | C      |  74.551 | 63.7204 | 57959.9469329 | 811.366817538894
      29 | C      | 67.9591 | 62.3363 | 57959.9553125 | 747.801892040569
      30 | C      | 68.9761 | 64.2901 | 57959.9639236 | 730.588787651341
     113 | C      | 86.0136 | 58.6753 | 57973.8089815 | 1038.88639325479
    (7 rows)

- Get all matched catalog sources by first getting the names of the
  catalogs and then querying for the source id number in those catalogs::

    SELECT * FROM catalog_match WHERE assoc_id = 1482;
    
      id  | catalog_id | src_id | assoc_id | separation 
    ------+------------+--------+----------+------------
     3275 |         12 | 499726 |     1482 |    8.16791
     3293 |         13 | 919076 |     1482 |    8.12618
    (2 rows)

    
    SELECT * FROM radcat.catalogs WHERE id IN (12, 13);
    
     id |   name   | telescope | frequency | resolution |      reference       
    ----+----------+-----------+-----------+------------+----------------------
     12 | nrl_nvss | VLA       |      1400 |         45 | 
     13 | nvss     | VLA       |      1400 |         45 | Condon et al. (1998)
    (2 rows)

    SELECT * FROM radcat.nrl_nvss WHERE id = (
        SELECT src_id FROM catalog_match WHERE assoc_id = 1482 AND
	    catalog_id = 12);
    
    SELECT * FROM radcat.nvss WHERE id = (
        SELECT src_id FROM catalog_match WHERE assoc_id = 1482 AND
	    catalog_id = 13);

- Compare images which contain a VLITE unique source, detected or not::

    SELECT image.id, filename, obs_ra, obs_dec, bmaj, bmin, noise,
        config, nvis, tau_time, nsrc, error_id, nearest_problem, separation,
        assoc_id, detected
    FROM image JOIN vlite_unique vu ON image.id = vu.image_id
    WHERE assoc_id = 10966 ORDER BY image.id;


.. _admin:

Database Admin Tasks
--------------------
Try not to break anything.

Changing the Data Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. highlight:: none
		
Follow these steps if the database cluster needs to be moved to a
different file location in the future. This has been tested once
without incident, so it should hopefully work for you. Execute the
commands as the postgres superuser unless otherwise specified.

1. Confirm the location of the current data directory. This
   should be whatever the PGDATA environment variable is set to
   for the vpipe and postgres users on roadrunner::

     $ echo $PGDATA
     $ /usr/local/pgsql/data
     
2. Shutdown the database server::

     $ pg_ctl stop
     
3. Make sure the database server is actually shutdown. *You will corrupt
   all data and ruin everything if it isn't*::

     $ pg_ctl status
     
4. Create the new data parent directory and give the
   postgres user ownership (*execute as vpipe*)::

     $ mkdir /data3/vpipe/vdp
     $ sudo chown postgres /data3/vpipe/vdp
     
5. Log back in as the postgres superuser. Copy the existing
   data directory to the new location using rsync::

     $ rsync -av /usr/local/pgsql/data /data3/vpipe/vdp

   This will create the directory /extraid/vpipe/vdp/data with the
   original permissions and postgres user ownership of the data
   directory.
   
6. Rename the old data directory to avoid any potential PostgreSQL
   confusion::

     $ mv /usr/local/pgsql/data /usr/local/pgsql/data_old

7. Point PGDATA to the new location. Update the variable for both
   vpipe and postgres users.
   
8. Start up the database server and verify that all the data
   is there and everything still works as expected. You can
   verify that PostgreSQL is accessing the correct data directory
   by executing these SQL statements as superuser postgres
   after the server has been started::

     $ psql
     postgres=> SHOW data_directory;
     postgres=> \q

9. You can remove the old data directory once you're confident
   that everything is working correctly with the new one::

     $ rm -rf /usr/local/pgsql/data_old

   *Be careful not to delete the pgsql parent directory.*
   All the executable files are in there, so you would probably
   have to re-install PostgreSQL.

Managing Users
^^^^^^^^^^^^^^
Database roles have already been created for most people
and privileges to the VLITE databases have been
set such that everyone has read-only access (SELECT statements
only). This even includes any new roles that might get added later.
Below are the steps that were taken to make this happen. All
steps must be completed as the postgres superuser.

1. Create a new database user/role::

     $ su - postgres
     $ createuser --interactive [username]

   [username] should be the same as the name of the person's
   account on roadrunner (i.e. erichards). Select no for all
   permissions.

2. Connect to each VLITE database and execute the following
   SQL statements to give every role read access to the database
   tables:

   .. code-block:: sql

	GRANT SELECT ON ALL TABLES IN SCHEMA public TO PUBLIC;
	GRANT USAGE ON SCHEMA radcat TO PUBLIC;
	GRANT SELECT ON ALL TABLES IN SCHEMA radcat TO PUBLIC;

   PUBLIC is a special role name that is used to grant a privilege to
   every role on the system.

Adding quantile extension
^^^^^^^^^^^^^^^^^^^^^^^^^
Tomas Vondra's quantile extension is useful and can be added to a database by running::

     $ psql dbname -c "CREATE EXTENSION quantile"

