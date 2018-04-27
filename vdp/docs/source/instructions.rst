.. _instructions:

Maintainer Instructions
=======================
This page contains how-to guides for various updates and
improvements to **vdp** which may have to be done in the future.

.. _add_new_catalog:

Adding a New Radio Catalog
--------------------------
Adding a new catalog of radio sources to be used for cross-matching
against VLITE sources will require editing the ``radiocatalogs.catalogio``
and ``radiocatalogs.radcatdb`` modules. Properties of each catalog
source are read in and stored as attributes of objects from the
``radiocatalogs.catalogio.CatalogSource()`` class. These objects and their
attributes are written to a text file which is then read and bulk
inserted into the database as a new table in the "radcat" schema.
Follow the instructions below to add or update a catalog.

1. If adding a new catalog, create a new entry for it in the
   catalog dictionary at the top of the ``radiocatalogs.catalogio``
   module. Add the name of the catalog as a new key and a
   dictionary as its value which contains the key 'id' and a
   new id number value. The catalog name must be in all lowercase,
   cannot lead with a number, and cannot contain the "." character.
   If you are updating and replacing an existing catalog, then you
   do not need to do anything to the catalog dictionary, unless
   you want to update its name.

   .. note:: Add the new catalog to the end of the dictionary with
	     the correct corresponding id number. If you insert it
	     elsewhere (i.e. to maintain alphabetical order) and
	     re-assign the id values of the other existing catalogs,
	     then you will need to re-create the "radcat" schema
	     *AND* re-write the "\*_psql.txt" files so that the
	     'catalog_id' of the sources matches the 'id' of the
	     catalog in the **radcat.catalogs** table.

2. Write a new function, or update the existing, to read the
   catalog in the ``radiocatalogs.catalogio`` module. Name the
   function "read_[catalog name]", where "[catalog_name]" is
   whatever you specified in the catalog dictionary. Use one
   of the other catalog reading functions as a template.
3. The final step is to add a couple of lines, or update the
   existing ones if necessary, in the ``radiocatalogs.radcatdb.add_table``
   function to call the function to read the catalog::

     elif tblname == '[catalog name]':
         sources = catalogio.read_[catalog name]()

.. _update_qa:

Updating the Image & Source Count Quality Checks
------------------------------------------------
The default values for the quality check requirements can
be updated in the ``vdp.cfgparse`` function. The code which
performs the checks is located in the ``database.dblcasses``
module under the ``Image`` class methods ``image_qa`` and
``source_qa``. If a new quality check is added, the
``database.createdb.make_error`` function will need to
be updated to add the new requirement to the database
**error** table. Additional parameters may also be added to
the configuration file along with code to parse them and set
default values in ``vdp.cfgparse``.

.. _maintenance:

Database Maintenance
--------------------
The ``maintenance.py`` code should be run periodically (frequency
TBD) once the database is large. Like ``vdp.py``, ``maintenance.py``
is called from the command line and will connect to the database
specified in the configuration file. It then runs 'CLUSTER' and
'ANALYZE' SQL commands on the **detected_source** and **assoc_source**
tables. Clustering orders the data on the disk according to the Q3C
spatial index values and analyze collects statistics about the contents,
which is used by the PostgreSQL query planner. Together, these procedures
ensure faster queries.
