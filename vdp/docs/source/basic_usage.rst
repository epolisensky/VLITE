.. _basic_usage:

Basic Usage
===========

**vdp** requires a YAML configuration file which tells it
where to find the VLITE images and which local database to
connect to (see :ref:`config_desc` for details).
Since it was built specifically for VLITE, **vdp** assumes
a particular directory structure in which it looks for all
files ending with "IPln1.fits":

  ``/root_directory/year-month/day/Images/``

Or, for example: ``/extraid/vpipe/processed/2018-03/26/Images/``

If you leave the "month" or "year" parameters blank in the
configuration file, **vdp** will look for VLITE images
in ``/root_directory/Images/``.

.. note:: The user will need write permissions in ``Images/``
	  and its parent directory.

The database can be created at runtime, so it does not need
to exist before starting the pipeline. If connecting to an
existing database, the user must either be the owner of that
database or a superuser due to the permissions needed while
running **vdp**. Unless you have a specific reason not to,
it is best to connect as user 'vpipe', which has already
been created on the PostgreSQL server running on virgo.

Running the Code
^^^^^^^^^^^^^^^^
Once the YAML configuration file is set, **vdp** can
be started from the command line::
  
  $ python vdp.py config_file.yaml

where ``config_file.yaml`` is replaced with whatever
you named your configuration file. The pipeline's progress
will be displayed through messages printed to the terminal.

*******************************
Optional Command Line Arguments
*******************************
There are some optional command line arguments that enable
some non-standard functionality for **vdp**.
Use the ``--help`` or ``-h`` command line option to see all
optional command line arguments::
  
  $ python vdp.py -h
  usage: vdp.py [-h] [--ignore_prompt] [--remove_catalog_matches]
           [--remove_source] [--add_catalog]
           config_file

  Run the VLITE Database Pipline (vdp)

  positional arguments:
    config_file           the YAML configuration file

  optional arguments:
    -h, --help            show this help message and exit
    --ignore_prompt       ignore prompt to verify database removal/creation
    --remove_catalog_matches
                          remove matching results for the specified sky survey
                          catalog(s)
    --remove_source       removes the specified source(s) from the database
                          assoc_source table
    --add_catalog         adds any new sky survey catalogs to a table in the
                          database "skycat" schema

``--ignore_prompt`` overrides the prompt to confirm overwriting
or creating a new database. This is needed when redirecting all
terminal output to a log file, as is recommended when batch
processing (see :ref:`batch_proc`).

``--remove_catalog_matches`` will prompt the user to input the
name(s) of all radio catalogs whose matching results are to be
removed from the databse. This will likely only be used in a
scenario where an updated version of a radio catalog is released
and you wish to replace the old matching results with new ones.

``--remove_source`` prompts the user for a list or text file
containing a list of id numbers corresponding to the rows to
remove in the database **assoc_source** table. Use this when
you want to remove sources you have determined to be artifacts
from the VLITE catalog. These sources will remain in the
database **detected_source** table, but are given an 'assoc_id'
value of -1.

``--add_catalog`` cross-checks the list of available radio
catalogs defined in ``skycatalog.catalogio.catalog_list``
with the table names that exist in the database "skycat"
schema. If a catalog is found in the list that does not
yet exist in the "skycat" schema, a new table is created
for it and the **skycat.catalogs** table is updated accordingly.
It will only be necessary to run this option after you have
added code to ``skycatalog.catalogio`` and
``skycatalog.skycatdb`` to read a new radio catalog and need
to add it to an existing database's "skycat" schema.
See :ref:`add_new_catalog` for more information.

.. _batch_proc:

****************
Batch Processing
****************
The configuration file enables processing of one ``year-mo``
directory at a time.
Processing more than one month of VLITE images can be accomplished
through successive runs of **vdp** called from a bash script.
It is recommended that all output that normally gets printed
to the terminal window be redirected to a text file to retain
a record of the pipeline's execution. Don't forget to use the
optional command line argument ``--ignore_prompt`` for the
first call to **vdp** if creating a new database or overwriting
an existing one.

Example file ``batch_vdp.bash``:
::
   
  python vdp.py 201801config.yaml --ignore_prompt > 201801.log
  python vdp.py 201802config.yaml > 201802.log
  python vdp.py 201803config.yaml > 201803.log

************************
Expected Execution Times
************************
Execution time mostly depends on the number and size of the
images being processed. Expect ~30-60 seconds per image for
VLA A configuration, 15-45 s/image for B, and 5-15 s/image
for C & D configurations, on average. The bottleneck is source
finding/measurement with PyBDSF.

*************
Data Products
*************
A ``PyBDSF/`` directory is created in the ``Images/`` parent directory
which stores the PyBDSF generated log files and ds9 regions
files for each image. The database contains all results from
each stage of the pipeline. See :ref:`database` for more
information.


.. _config_desc:

Description of Configuration File Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example of the required YAML configuration file can be
found in the VLITE GitHub repository `here.
<https://github.com/erichards/VLITE/blob/develop/vdp/example_config.yaml>`_
The contents are described in more detail below.

**stages**
  Accepts boolean ``True``/``False`` or "yes"/"no" to turn on/off
  running certain pipeline stages.

  *source finding*
    Runs source finding & measurement on the image with PyBDSF.
    (See :ref:`source_finding`).
  *source association*
    Associates the image's detected sources with the existing VLITE
    catalog contained in the database **assoc_source** table.
    (See :ref:`source_assoc`).
  *catalog matching*
    Cross-matches the image's detected sources with sources from
    other radio catalogs.
    (See :ref:`catalog_matching`).

**options**
  Accepts boolean ``True``/``False`` or "yes"/"no" to turn on/off
  certain features for the pipeline.

  *save to database*
    Saves all results to the database.
  *quality checks*
    Checks if the image meets certain quality standards before
    and after source finding. (See :ref:`image_qa` and
    :ref:`source_count_qa`).
  *overwrite*
    Deletes all contents & re-creates tables, functions, triggers,
    and indices in the existing database "public" schema.
  *reprocess*
    Deletes all existing results for the image and re-runs source
    finding plus any additional stages specified. Applies only
    if the source finding stage is turned on.
  *redo match*
    Deletes all matching results between the image's detected
    sources and other radio catalogs. Cross-matching is then
    run again for those image's sources using the currently
    specified list of radio catalogs.
  *update match*
    Cross-matches the image's detected sources with any currently
    specified radio catalogs for which there are no results yet.
  
**setup**
  Parameters defining location of VLITE images and database
  connection info.

  *root directory*
    Root path to the VLITE images (i.e. ``/extraid/vpipe/processed/``).
  *year*
    Four-digit calendar year (i.e. ``2018``). If blank, directory
    path is ``/root_directory/Images/``
  *month*
    One- or two-digit numerical calendar month (i.e. ``03``).
    If blank, directory path is ``/root_directory/Images/``
  *day*
    List of two-digit daily directories to process under the
    ``year-mo`` parent directory. To process all, leave as
    empty list, ``[]``. Otherwise, ``[01, 02, 03, etc.]``.  
  *database name*
    Name of new or existing database.
  *database user*
    Name of the PostgreSQL database user.
  *catalogs*
    List of other radio catalogs to use for cross-matching. To use all
    available catalogs, leave as empty list, ``[]``.
    Otherwise, ``[FIRST, TGSS, NVSS, WENSS, VLSSr, etc.]``.

**pybdsf_params**
  Parameters used in source finding.

  *mode*
    Required -- choose either 'default' or 'minimize_islands'.
    Determines whether PyBDSF is run once per image ('default'),
    or multiple times with different ``rms_box`` parameters to
    find the fewest number of islands ('minimize_islands').  
  *scale*
    Required -- number between 0 and 1. Fraction of the image's
    field-of-view to use. The length of the radius describing
    the image's circular field-of-view is multiplied by this number.

  Below this point, any number of PyBDSF parameters may be
  specified. See `their documentation <http://www.astron.nl/citt/pybdsm/process_image.html#general-reduction-parameters>`_ for descriptions of
  all available options. The parameters shown below have been
  found to work best for VLITE images:
  
    - ``thresh``: 'hard'
    - ``adaptive_rms_box``: ``True``
    - ``adaptive_thresh``: 10.
      
  If you want to specify any PyBDSF parameter that accepts a
  tuple, like ``rms_box``, it needs to be formatted as such:
  
    rms_box: !!python/tuple [100, 30]
  
**image_qa_params**
  Sets quality requirements for images. Applies only if quality checks
  are turned on. Leave any parameter blank to use the default value.

  *min time on source (s)*
    Minimum allowed integration time on source. Image header
    keyword ``TAU_TIME``. Default is 60 seconds.
  *max noise (mJy/beam)*
    Maximum allowed image noise. Image header keyword ``ACTNOISE``.
    Default is 1000 mJy/beam.
  *max beam axis ratio*
    Maximum allowed ratio between the beam semi-major and
    semi-minor axes. Default is 4.
  *min problem source separation (deg)*
    Minimum allowed angular separation between the image
    pointing center and a known problem source/area.
    Default is 20 degrees.
  *max source metric*
    Maximum allowed metric for source counts. Defined as:
    (actual_num_sources - expected_num_sources) / expected_num_sources.
    Default is 10.
