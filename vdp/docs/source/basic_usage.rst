.. _basic_usage:

Basic Usage
===========

**vdp** requires a YAML configuration file which tells it
where to find the VLITE images and which local database to
connect to (see :ref:`config_desc` for details).
Since it was built specifically for VLITE, **vdp** assumes
a particular directory structure in which it looks for all
files ending with "IPln1.fits" for VLITE or "IMSC.fits" for VCSS:

  ``/root_directory/year-month/day/image_directory/``

Or, for example, ``/extraid/vpipe/processed/2018-03/26/Images/``
for VLITE or, ``/data3/vpipe/vcss/2017-09/21/Mosaics/`` for VCSS.

If you leave the "month" or "year" parameters blank in the
configuration file, **vdp** will look for VLITE images
in ``/root_directory/Images/``.

.. note:: The user will need write permissions in image directory
	  and its parent directory.

The database can be created at runtime, so it does not need
to exist before starting the pipeline. If connecting to an
existing database, the user must either be the owner of that
database or a superuser due to the permissions needed while
running **vdp**. The database user vpipe should be used for
running **vdp** on roadrunner.

Running the Code
^^^^^^^^^^^^^^^^
Once the YAML configuration file is set, **vdp** can
be started from the command line::
  
  $ python vdp.py config_file.yaml

where ``config_file.yaml`` is replaced with whatever
you named your configuration file. The pipeline's progress
will be displayed through messages printed to the console
and a log file.

**********************
Command Line Arguments
**********************
There are some optional command line arguments that enable
some non-standard functionality for **vdp**.
Use the ``--help`` or ``-h`` command line option to see all
required and optional command line arguments::
  
  $ python vdp.py -h
  usage: vdp.py [-h] [--ignore_prompt] [--remove_catalog_matches]
           [--remove_source] [--add_catalog]
           config_file

  Run the VLITE Database Pipline (vdp)

  positional arguments:
    config_file           the YAML configuration file

  optional arguments:
    -h, --help            show this help message and exit
    -q, --quiet           stops printing of messages to the console
    --ignore_prompt       ignore prompt to verify database removal/creation
    --remove_catalog_matches
                          remove matching results for the specified sky survey
                          catalog(s)
    --remove_source       removes the specified source(s) from the database
                          assoc_source table
    --remove_image        removes the specified image(s) and associated results
                          from the database entirely
    --manually_add_match  manually add catalog matching results for VLITE
                          source(s) after follow-up
    --add_catalog         adds any new sky survey catalogs to a table in the
                          database "skycat" schema
    --update_pbcor        update corrected_flux table with primary beam
                          corrections

``-q, --quiet`` turns off printing of log messages to
the console. They will still be recorded in the log file.

``--ignore_prompt`` overrides the prompt to confirm overwriting
or creating a new database. This may be desired when batch
processing (see :ref:`batch_proc`) and suppressing console
print statements with ``-q``.

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

``--remove_image`` prompts the user for a list or text file
containing the filenames of the images to delete from the database.
All results from the specified images are removed from every
affected table. The image filenames must contain the full
directory path starting at least from the year-month directory
structure.

``--manually_add_match`` enables the user to add matching results
between VLITE and other radio catalog sources that failed to be
matched due to the angular separation being just over the limit.
The id of the VLITE source in the **assoc_source** table and the
name of the catalog must be provided at a minimum. The catalog
source id and angular separation may also be included. The required
inputs can be provided on the command line, one assoc_id, catalog
name combination at a time, pressing "q" when finished.
Alternatively, a text file containing the necessary information
may be given.

``--add_catalog`` cross-checks the list of available radio
catalogs defined in ``radiocatalogs.catalogio.catalog_list``
with the table names that exist in the database "radcat"
schema. If a catalog is found in the list that does not
yet exist in the "radcat" schema, a new table is created
for it and the **radcat.catalogs** table is updated accordingly.
It will only be necessary to run this option after you have
added code to ``radiocatalogs.catalogio`` and
``radiocatalogs.radcatdb`` to read a new radio catalog and need
to add it to an existing database's "radcat" schema.
See :ref:`add_new_catalog` for more information.

``--update_pbcor`` reads each row of the image table, fetches
each image's detected_sources and updates the corrected_flux table
with primary beam corrections. Intended for use when new primary
beam models are available. Beam models for each primary observing
frequency are set in read_fitted_beam in sourcefinding/beam_tools.py

.. _batch_proc:

****************
Batch Processing
****************
The configuration file enables processing of one ``year-mo``
directory at a time.
Processing more than one month of VLITE images can be accomplished
through successive runs of **vdp** called from a bash script.
You can suppress output to the console by using ``-q`` or
``--quiet``. All output will be written to a log file
in the root directory with name "yearmo.log" (i.e. "201801.log").
You may additionally use the optional command line argument
``--ignore_prompt`` for the first call to **vdp** if creating
a new database or overwriting an existing one and don't want to
stick around to verify.

Example file ``batch_vdp.bash``:
::
   
  python vdp.py 201801config.yaml -q --ignore_prompt
  python vdp.py 201802config.yaml -q
  python vdp.py 201803config.yaml -q

************************
Expected Execution Times
************************
Execution time mostly depends on the number and size of the
images being processed. Expect ~45-90 seconds per image for
VLA A configuration, 15-45 s/image for B, and 5-15 s/image
for C & D configurations, on average. The bottleneck is source
finding/measurement with PyBDSF.

*********************
Stopping the Pipeline
*********************
Execution times can be long (hours/days) when processing many
images. There may be times when you need to stop the pipeline
before it has completed and restart it later. A handler has
been implemented (thanks to an internet blogger) to gracefully
break out of the processing loop. A keyboard interrupt (CTRL-C)
will signal the pipeline to stop once it has finished processing
the current image and exit as if the run had completed normally.

You can simply restart the pipeline with same configuration file.
**vdp** will skip any file it finds is already in the database
*image* table if the *reprocess* configuration file option is
turned off. You may also want to edit the *day* and/or *files*
lists in the configuration file to run only the ones remaining
so there aren't hundreds of lines printed about skipping over
already-processed images.

If things have gone completely off the rails and you need
to kill the pipeline NOW, hitting CTRL-C nine times will
override the graceful exit and send a real keyboard interrupt
to Python. Basically, just keep doing CTRL-C until everything
comes grinding to a halt. No guarantees on the state of the
database after that, though.

*************
Data Products
*************
A ``PyBDSF/`` directory is created in the image parent directory
which stores the PyBDSF generated log files and ds9 regions
files for each image. A log file is also generated in the root
directory, or appended to if it already exists, with every run of
the pipeline. The database contains all results from
each stage of the pipeline. See :ref:`database` for more
information.


.. _config_desc:

Description of Configuration File Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example of the required YAML configuration file can be
found in the VLITE GitHub repository `here.
<https://github.com/erichards/VLITE/blob/master/vdp/example_config.yaml>`_
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
  *beam corrected*
    Are the images already primary beam corrected?
  *always associate*
    Associates sources in all images regardless of image ass_flag
  
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
  *image directory*
    Name of the sub-directory where the image files are
    located under ``root_directory/year-month/day``.
    The default is ``Images/`` if left blank.
  *files*
    Lists of files to process in each daily directory. To process
    all, leave as empty nested list, ``[[]]``. Otherwise,
    ``[[f1.fits, f2.fits, etc.], [f1.fits, etc.], etc.]``
  *database name*
    Name of new or existing database.
  *database user*
    Name of the PostgreSQL database user.
  *catalogs*
    List of other radio catalogs to use for cross-matching. To use all
    available catalogs, leave as empty list, ``[]``.
    Otherwise, ``[FIRST, TGSS, NVSS, WENSS, VLSSr, etc.]``.
  *smeartime*
    Max time step between primary beam sampling in each continuous 
    time interval (from NX table in UVOUT file). (default: 900 [s])

**pybdsf_params**
  Parameters used in source finding.

  *mode*
    Required -- choose either 'default' or 'minimize_islands'.
    Determines whether PyBDSF is run once per image
    ('default'; recommended), or multiple times with different
    ``rms_box`` parameters to find the fewest number of islands
    ('minimize_islands').  
  *scale*
    Required -- number between 0 and 1. Fraction of the image's
    field-of-view to use. The length of the radius describing
    the image's circular field-of-view is multiplied by this number.
  *borderpad*
    Required -- sources within this many pixels of the image border 
    will be rejected. Default value is 3

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

  *min nvis*
    Minimum allowed number of visibilities. Image header
    keyword ``NVIS``. Default is 1000 seconds.
  *max sensitivity metric*
    Maximum allowed combination of noise & integration time on source.
    Defined as noise x sqrt(int. time). Default is 3000 mJy/beam s^1/2.
  *max beam axis ratio*
    Maximum allowed ratio between the beam semi-major and
    semi-minor axes. Default is 4.
  *max source count metric*
    Maximum allowed metric for source counts. Defined as:
    (actual_num_sources - expected_num_sources) / expected_num_sources.
    Default is 10.
  *min niter*
    Minimum number of CLEAN interations. Important for 
    reliable source fluxes. Image header keyword ``NITER``
    or ``CLEANNIT`` or ``NITER`` in a ``HISTORY`` line. 
    Default is 1000.
  *min bpix*
    Minimum size of BMIN in pixels (defualt: 2.8)
  *max bpix*
    Maximum size of BMIN in pixels (default: 7)
