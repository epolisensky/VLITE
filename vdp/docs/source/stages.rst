.. _stages:

Pipeline Stages
===============
Data flow through the VLITE Database Pipeline can be broken into four
stages. These stages can be run in succession on one image at a time
in a single execution of the pipeline, or they can be run one stage
at a time in multiple executions of the pipeline while processing
multiple images each time in each separate stage. The former is the
default preferred method since the latter requires multiple executions
of the pipeline. Stages can be turned on and off through the YAML
configuration file. They are described in more detail below.

.. _read_image:

Stage 1: Reading the Image
--------------------------
The first stage of the pipeline simply reads the VLITE image and
records some of its header metadata in the database **image** table.
This stage is not directly specified in the configuration file, but
is always executed first since it is necessary for all other stages.
Images are processed in time order based on their Modified Julian
Date.

**vdp** will skip processing an image if it sees that its filename
(including the full directory path) is already in the **image** table
and the image has already been run through the stages specified in the
configuration file. If you wish to force re-processing of an image
that's already in the database, use the configuration file option
*reprocess*.

.. note:: The *reprocess* option removes all previous source finding
	  results from the database and only applies if the
	  source finding stage is turned on.

If the *save to database* option is turned off, then the image header
info will not be inserted into the **image** table. The image is
simply read into Python to be passed along to the other stages. If
all stages are additionally turned off, **vdp** will exit with a
``ConfigError`` complaining that there's nothing to do.

.. _image_qa:

Image Quality Assurance
^^^^^^^^^^^^^^^^^^^^^^^
VLITE data is vastly inhomogeneous, so quality assurance is
important to filter out unusable images. The first round
of quality checks attempts to flag bad images based on their
header metadata. Images are checked for the following requirements
in the order they are listed:

1. the image needs to contain all header keywords used in the
   pipeline, with the exception of ``OBJECT``, ``DATE-MAP``,
   ``PEAK`` or ``DATAMAX``, ``CONFIG``, and ``DURATION``
2. the integration time on source (``TAU_TIME``) must be
   greater than some minimum time (default is 60 s)
3. the image noise must be >= 0 and <= some maximum value
   (default is 500 mJy/beam)
4. the ratio of the image's beam semi-major to semi-minor axis
   must be less than some maximum value (default is 4)
5. the imaging target can't be a planet or the NCP
6. a known problem source cannot be in the image's field-of-view:
   - the Sun
   - Cassiopeia A
   - Cygnus A
   - Taurus A
   - Hercules A
   - M87 (Virgo A)
   - Galactic Center

Each of the above requirements has a default value which is used
unless otherwise specified through the corresponding option under
**image_qa_params** in the configuration file. The above requirements
and their values are recorded in the database **error** table.

.. note:: The **error** table records the values specified under
	  **image_qa_params**, or the defaults if left blank, in
	  the configuration file *at the time
	  of creation*. The table is not updated if the requirements
	  change in subsequent runs using the same database. They
	  are, however, recorded in the **run_config** table and
	  output at runtime if an image fails a quality check.

If an image fails a quality check, it is assigned an *error_id*
corresponding to the requirement it failed to pass which can be
looked up in the **error** table. Once an image fails a check,
its information is added to the **image** table and the pipeline
moves on to the next image. It does not push the failed image
through the remaining quality checks, if any, or allow it to pass
through to any other stages. Quality checks can be turned off
in the configuration file.

.. _source_finding:

Stage 2: Source Finding
-----------------------
**vdp** integrates PyBDSF to automate finding and measuring
sources in the VLITE images. PyBDSF measures sources by grouping
Gaussians that were fit to islands of contiguous pixels.
Islands are formed by finding all pixels
with a peak flux greater than 5-sigma and includes all surrounding
pixels above 3-sigma. These threshold values represent the defaults
and can be changed under **pybdsf_params** in the configuration
file. Identification of pixels above the thresholds depends on the
local mean and rms. PyBDSF calculates background mean and rms images
using a sliding 2-D box with interpolation. The size of this box
and the step size inbetween are controlled by the ``rms_box`` parameter.
Refer to the PyBDSF `documentation
<http://www.astron.nl/citt/pybdsm/index.html>`_ for more information.

Properties of the sources and islands are written to the database
**detected_source** and **detected_island** tables, respectively.
A ds9 region file is also created for every image. A 1-D primary
beam correction factor is applied to all flux measurements from
PyBDSF and recorded in the **corrected_flux** table. A 20%
systematic uncertainty is also added in quadrature to the PyBDSF
reported 1-sigma statistical uncertainties on all flux error
measurements in the **corrected_flux** table only. The applied
primary beam correction factor was determined empirically and
depends only on the source's distance from the image center,
which is also recorded in the **corrected_flux** table in degrees.
Corrected fluxes are only computed if the *save to database*
option is turned on in the configuration file.

**Pre-Processing**

The VLITE images are cropped to a circular field-of-view before being
ingested by PyBDSF. This is done to ensure that cone search queries
in the database return sources which lie in the same well-defined
field-of-view as the images. The radius of the circular field is half
the image size, which for VLITE is 1 degree for A cofiguration, 2
degrees for B/B+, 3 degrees for C, and 4 degrees for D. The *scale*
parameter in the **vdp** configuration file can be used to make the
field-of-view radius smaller by the factor specified.

Experience has shown that it is better to specify the ``rms_box``
parameter for VLITE images rather than leave PyBDSF to calculate it
internally (the default option). Images with bright stripe artifacts
will fail miserably with the default option. Setting the box size to
1/10th of the image size in pixels and the step size to 1/3rd of the
box size seems to help avoid identifying bright non-point-like
artifacts as real emission while still capturing most of the real
point sources. The ``rms_box`` parameter is calculated in this way
for every image run through **vdp**, unless it is explicitly specified
in the configuration file. The ``rms_box_bright`` parameter is also
calculated for every image as 1/5th times the rms_box sizes. This
parameter is used by PyBDSF only if ``adaptive_rms_box`` is set to
``True`` in the configuration file. This tells PyBDSF to use the
smaller box size around bright sources where there tend to be more
imaging artifacts. See the :ref:`config_desc` for additional PyBDSF
parameters to specify explicitly for VLITE under **pybdsf_params**
in the configuration file.

.. _source_count_qa:

Source Count Quality Assurance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A second round of quality checks are performed on the source finding
results (again, only if the *quality checks* option is turned on)
before they are inserted into the database tables. Images are flagged
if PyBDSF failed to process it for any reason or if there were no
sources extracted. Any image that takes longer than 5 minutes for
PyBDSF to process will fail with a timeout error to avoid PyBDSF
getting stuck trying to fit Gaussians to large imaging artifacts.
We also define a metric developed by E. Polisensky to flag images
where the number of detected sources is much larger than what is
expected based on source counts from the WENSS survey and the image's
noise. The absolute difference between the actual number of sources and the
expected number of sources normalized by the expected number is
required to be less than some value (default is 10).

As with the initial image quality checks, images that fail
will be assigned an *error_id* corresponding to a PyBDSF failure to
process, zero sources found, or an unrealistic number of sources
found which is recorded in the **image** table. These image's
sources, if there are any, are not carried forward to the
association or catalog matching stage and are not written to
the **corrected_flux** table.

.. _source_assoc:

Stage 3: Source Association
---------------------------
The association stage condenses multiple detections of a single source
from different images into one entry in the **assoc_source** database
table. Detections of the same source are required to be at similar
spatial resolutions before being associated to avoid differences in
source structure (i.e. resolved double vs. unresolved single). The
resolution of an image is defined by the beam semi-minor axis size so
it is less sensitive to elongated beam shapes. Currently, images are
divided into four resolution classes which roughly correspond to the
four VLA configurations:

- resolution <= 15" (A/VLITE B+)
- 15" < resolution <= 35" (B)
- 35" < resolution <= 60" (C)
- 60" < resolution (D)

After source finding, the association stage proceeds as follows:

1. A cone search query is sent to the **assoc_source** table to extract
   all sources detected in previous images which lie in the same
   field-of-view as was used in source finding on the current image.
2. The extracted association candidates are then filtered on 'res_class'
   so that only candidates in the same resolution class as the current
   image remain.
3. Sources detected in the current image are cross-matched with the
   filtered association candidates by choosing nearest neighbors that are
   separated by less than half the length of the current image's beam
   semi-minor axis.
4. If a successful association is made, the position of the source
   recorded in the **assoc_source** table is updated to reflect the
   weighted average of all detections and the number of detections,
   'ndetect', is incremented. If no association is made, those detected
   sources are added to the **assoc_source** table as new sources in
   that resolution class.
5. A new entry in the **vlite_unique** table is added for every
   association candidate pulled from the **assoc_source** table with no
   catalog matches ('nmatches' = 0) to record another VLITE detection
   for that source in the current image if there was a successful
   association or to record a non-detection for the current image if
   there was not.

If the *save to database* option is turned off, the association results
are printed to the console and/or log file without updating the
database tables.

.. _catalog_matching:

Stage 4: Catalog Matching
-------------------------
All VLITE sources are cross-matched with other radio sky surveys and
catalogs to help isolate transient candidates and compare fluxes across
the radio spectrum. As for the association stage, cross-matching is
restricted between sources with similar spatial resolutions -- the
resolution of the catalog has to be in the same resolution class as the
image. The resolution classes are the same as for association except the
first two classes (A & B) are combined.

The cross-matching steps proceed as follows:

1. The list of catalogs specified in the configuration file is filtered
   to remove ones outside the acceptable range of spatial resolution for
   the current image.
2. The 'catalogs_checked' column in the **image** table is queried to
   see which, if any, of the resolution-filtered catalogs have already
   been checked for the current image. Only new catalogs which have not
   yet been checked for matches are used going forward.
3. VLITE sources are cross-matched with sources from each new,
   resolution-filtered catalog using the same method as for association:
   nearest neighbors with a separation less than half the beam's
   semi-minor axis length.
4. If a match is successful, the id of the VLITE source in the
   **assoc_source** table is added to the **catalog_match** table along
   with the matched catalog source's id and catalog id. The number of
   catalog matches, 'nmatches', in the **assoc_source** table is
   incremented for the matched VLITE source. If no match is found,
   'nmatches' is set to 0 and the **assoc_source** id & **image** id of
   the VLITE source are added to the **vlite_unique** table. The
   **image** table is queried to find all previously processed images
   with the same field-of-view as the current image which would have
   contained the VLITE source. New entries are added to the
   **vlite_unique** table for those images to record their non-detection
   of the VLITE source.

Which VLITE sources are used in cross-matching depends on how the
pipeline is being run. When following the source association stage,
all new VLITE sources that are added to the **assoc_source** table
for the first time are passed on for catalog cross-matching. This is
so every associated VLITE source is only cross-matched with other radio
catalogs once. If the existing catalog matching results need to be
re-done for the current image, this can be accomplished by turning on
*redo match* in the configuration file. With *redo match* set to ``True``,
catalog matching results will be wiped clean for entries in the
**assoc_source** table that correspond to sources detected in the
current image and then re-matched with sources from the currently
specified list of catalogs. It is also possible to add a new catalog
to existing results without re-doing all cross-matching for all catalogs
by turning on *update match* in the configuration file. If *update match*
is ``True``, all entries in the **assoc_source** table that correspond
to sources detected in the current image will be cross-matched against
sources from any currently specified catalog for which there are not
already matching results for that VLITE associated source.

The functionality to directly cross-match all sources detected in the
current VLITE image with any specified radio catalogs, regardless of
spatial resolution, is enabled by running only the source finding and
catalog matching stages with *save to database* disabled. These results
will then be printed to the console, but are not saved to the database.
