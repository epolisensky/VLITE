.. _database:

Database Contents & Structure
=============================
The **vdp** database keeps a record of every VLITE image,
detected source, association between detections of the
same source, and match to sources from other radio catalogs.
The :ref:`tables` section below lists every table in the
database along with a description of the columns.
A spreadsheet view of the database organization can be
seen `here. <https://docs.google.com/spreadsheets/d/e/2PACX-1vR20qGzJ7U3hFBNYZ1IUJWcFpdlOmfjQKv_8pk6aRW7BuljZ6VGNWyagHnsMVkZ6_Y9-Dl1vEwNv8Bg/pubhtml>`_

.. _tables:

Tables
^^^^^^

.. _assoc_source:

************
assoc_source
************
This is the VLITE catalog. It contains a single entry for every
source detected by associating detections from multiple images
in the same spatial resolution class.

- *id*: unique row id number
- *ra*: weighted average of all right ascension measurements (degrees)
- *e_ra*: error on the right ascension weighted average (degrees)
- *dec*: weighted average of all declination measurements (degrees)
- *e_dec*: error on the declination weighted average (degrees)
- *res_class*: spatial resolution class (see :ref:`source_assoc`)
- *ndetect*: number of images in which the source was detected
- *ns*: number of PyBDSF code 'S' detections of the source (single component Gaussian and only source in island)
- *nc*: number of PyBDSF code 'C' detections of the source (single component Gaussian in island with other sources)
- *nm*: number of PyBDSF code 'M' detections of the source (multi-Gaussian source)
- *nmatches*: number of matches to other radio catalog sources
- *ave_total*: weighted average of corrected total flux light curve (mJy)
- *e_ave_total*: error on the weighted average of corrected total flux light curve (mJy) 
- *ave_peak*: weighted average of corrected peak flux light curve (mJy/beam)
- *e_ave_peak*: error on the weighted average of corrected peak flux light curve (mJy/beam)
- *v_total*: variability metric of corrected total flux light curve
- *v_peak*: variability metric of corrected peak flux light curve
- *eta_total*: variability significance metric of corrected total flux light curve
- *eta_peak*: variability significance metric of corrected peak flux light curve

.. _catalog_match:

*************
catalog_match
*************
This table keeps a record of each radio catalog source that
was successfully cross-matched with a VLITE source in the
**assoc_source** table.

- *id*: unique row id number
- *catalog_id*: references the *id* column of the **radcat.catalogs**
  table
- *src_id*: references the *id* column of the radio catalog
  table (i.e. **radcat.nvss**) which corresponds to the
  *catalog_id*
- *assoc_id*: references the *id* column of the **assoc_source**
  table; deletion of row in **assoc_source** table cascades to
  all **catalog_match** rows with that *assoc_id*
- *separation*: angular separation between the weighted average
  VLITE source position and the other radio catalog source
  position (arcseconds)

.. _corrected_flux:

**************
corrected_flux
**************
All flux measurements from PyBDSF are corrected for VLITE's primary
beam response. A scale factor is applied to account for
decreasing sensitivity with increasing distance from the image's
pointing center. A primary observing frequency specific systematic 
uncertainty is also added (in quadrature) to the PyBDSF statistical 
uncertainties for every flux measurement error, typically 3%.

- *src_id*: non-unique id number assigned by PyBDSF to each detected
  source in an image; uniquely referenced by combination of
  (*src_id*, *image_id*), same as in **detected_source** table
  (deletions cascade from that table)
- *isl_id*: non-unique id number assigned by PyBDSF to each detected
  island in an image; uniquely referenced by combination of
  (*isl_id*, *image_id*), same as in the **detected_source** and
  **detected_island** tables (deletions cascade from **detected_island**)
- *image_id*: unique number identifying the image; references the
  *id* column in the **image** table and the *image_id* columns
  in the **detected_source** and **detected_island** tables
- *total_flux*: primary beam corrected total integrated flux (mJy)
- *e_total_flux*: 1-sigma statistical + systematic error on the
  corrected total integrated flux (mJy)
- *peak_flux*: primary beam corrected peak flux density per beam (mJy/beam)
- *e_peak_flux*: 1-sigma statistical + systematic error on the
  corrected peak flux (mJy/beam)
- *isl_total_flux*: primary beam corrected total integrated flux density
  of the island in which the source is located (mJy)
- *isl_e_total_flux*: 1-sigma statistical + systematic error on the
  corrected island integrated flux (mJy)
- *isl_rms*: average background rms noise of the island (mJy/beam)
- *isl_mean*: avearage background mean value of the island (mJy/beam)
- *isl_resid_rms*: average residual background rms noise of the island
  (mJy/beam)
- *isl_resid_mean*: average residual background mean value of the
  island (mJy/beam)
- *distance_from_center*: angular separation between the source position
  and the image pointing center (degrees)
- *polar_angle*: angle east of north of source in image plane (degrees)
- *snr*: signal-to-noise ratio of the source detection; defined as
  (*peak_flux* - *isl_mean*) / *isl_rms*
- *compactness*: measures source spatial extent from its *snr* and flux ratio; C > 1 is point-like, < 1 is extended source
- *clean*: True if source was CLEANed
- *xpix*: Source x coordinate in image (pixels)
- *ypix*: Source y coordinate in image (pixels)
- *assoc_id*: references the *id* column of the **assoc_source**;
  rows with the same *assoc_id* value are associated detections
  of the same source in different images

.. _detected_island:

***************
detected_island
***************
PyBDSF finds contiguous groups of pixels with flux values above a certain
detection threshold which it identifies as islands. Sources are extracted
from these islands by fitting them with Gaussians. This table stores
all PyBDSF output pertaining to the islands detected in an image (see
`the PyBDSF documentation <http://www.astron.nl/citt/pybdsm/write_catalog.html#definition-of-output-columns>`_ for a complete list).

- *isl_id*: non-unique id number assigned by PyBDSF to each detected
  island in an image; uniquely referenced by combination of
  (*isl_id*, *image_id*) & deletions cascade to **corrected_flux** table
- *image_id*: unique number identifying the image; references the
  *id* column in the **image** table (deletions cascade from that table)
- *total_flux*: total integrated flux density of the island (mJy);
  PyBDSF output column *Isl_Total_flux*
- *e_total_flux*: 1-sigma statistical error on the total integrated
  flux density (mJy); PyBDSF output column *E_Isl_Total_flux*
- *rms*: average background rms noise of the island (mJy/beam);
  PyBDSF output column *Isl_rms*
- *mean*: average background mean value of the island (mJy/beam);
  PyBDSF output column *Isl_mean*
- *resid_rms*: average residual background rms noise of the island
  (mJy/beam); PyBDSF output column *Resid_Isl_rms*
- *resid_mean*: average residual background mean value of the island
  (mJy/beam); PyBDSF output column *Resid_Isl_mean*

.. _detected_source:

***************
detected_source
***************
Properties of sources detected in an image are derived from Gaussian
fits to islands of pixels. This table records PyBDSF output
pertaining to the sources formed from fitting Gaussians (see
`the PyBDSF documentation <http://www.astron.nl/citt/pybdsm/write_catalog.html#definition-of-output-columns>`_ for a complete list). Descriptions
are either adapted or straight from their documentation.

- *src_id*: non-unique id number assigned by PyBDSF to each detected
  source in an image; uniquely referenced by combination of
  (*src_id*, *image_id*) & deletions cascade to **corrected_flux** table
- *isl_id*: non-unique id number assigned by PyBDSF to each detected
  island in an image; uniquely referenced by combination of
  (*isl_id*, *image_id*) & deletions cascade to **corrected_flux** table
  & from **detected_island** table
- *image_id*: unique number identifying the image; references the
  *id* column in the **image** table
- *ra*: source right ascension (degrees); PyBDSF output column *RA*
- *e_ra*: error on the right ascension (degrees); PyBDSF output
  column *E_RA*
- *dec*: source declination (degrees); PyBDSF output column *DEC*
- *e_dec*: error on the declination (degrees); PyBDSF output column
  *E_DEC*
- *total_flux*: total integrated flux (mJy); PyBDSF output
  column *Total_flux*
- *e_total_flux*: 1-sigma statistical error on the total integrated
  flux; PyBDSF output column *E_Total_flux*
- *peak_flux*: peak flux density per beam of the source (mJy/beam);
  PyBDSF output column *Peak_flux*
- *e_peak_flux*: 1-sigma statistical error on the peak flux (mJy/beam);
  PyBDSF output column *E_Peak_flux*
- *ra_max*: right ascension of the maximum of the source (degrees);
  PyBDSF output column *RA_max*
- *e_ra_max*: 1-sigma statistical error on the right ascension of
  the maximum (degrees); PyBDSF output column *E_RA_max*
- *dec_max*: declination of the maximum of the source (degrees);
  PyBDSF output column *DEC_max*
- *e_dec_max*: 1-sigma statistical error on the declination of
  the maximum (degrees); PyBDSF output column *E_DEC_max*
- *maj*: the FWHM of the major axis of the source (arcsec);
  PyBDSF output column *Maj*
- *e_maj*: 1-sigma statistical error on the FWHM of the source
  major axis (arcsec); PyBDSF output column *E_Maj*
- *min*: the FWHM of the minor axis of the source (arcsec);
  PyBDSF output column *Min*
- *e_min*: 1-sigma statistical error on the FWHM of the source
  minor axis (arcsec); PyBDSF output column *E_Min*
- *pa*: position angle of the source major axis measured east
  of north (degrees); PyBDSF output column *PA*
- *e_pa*: 1-sigma statistical error on the source major axis
  position angle (degrees); PyBDSF output column *E_PA*
- *dc_maj*: the FWHM of the deconvolved major axis of the source
  (arcsec); PyBDSF output column *DC_Maj*
- *e_dc_maj*: 1-sigma statistical error on the FWHM of the source
  deconvolved major axis (arcsec); PyBDSF output column *E_DC_Maj*
- *dc_min*: the FWHM of the deconvolved minor axis of the source
  (arcsec); PyBDSF output column *DC_Min*
- *e_dc_min*: 1-sigma statistical error on the FWHM of the source
  deconvolved minor axis (arcsec); PyBDSF output column *E_DC_Min*
- *dc_pa*: position angle of the source deconvolved major axis
  measured east of north (degrees); PyBDSF output column *DC_PA*
- *e_dc_pa*: 1-sigma statistical error on the source deconvolved
  major axis position angle (degrees); PyBDSF output column *E_DC_PA*
- *code*: defines the source structure:
  
  - 'S' = a single-Gaussian source that is the only source in the island
  - 'C' = a single-Gaussian source in an island with other sources
  - 'M' = a multi-Gaussian source

- *assoc_id*: references the *id* column of the **assoc_source**;
  rows with the same *assoc_id* value are associated detections
  of the same source in different images


.. _detected_null

***************
detected_null
***************
COMING SOON. Properties of sources NOT detected in an image that should've been. 
Derived from forced-fitting at the location of each unmatched
associated source an island with a size approximately the imaging beam. 
A null is considered detected if the ratio of source average 
flux to fitted flux is > 5.

- *assoc_id*: references the *id* column of the **assoc_source**;
  rows with the same *assoc_id* value are associated null detections
  of the same source in different images
- *image_id*: unique number identifying the image; references the
  *id* column in the **image** table
- *ra*: source right ascension of the fit (degrees) 
- *dec*: source declination of the fit (degrees)
- *total_flux*: beam corrected total integrated flux, set from island *max_value* (mJy)
- *e_total_flux*: beam corrected 1-sigma statistical error on the total integrated
  flux, set from island *total_fluxE* (mJy)
- *peak_flux*: beam corrected peak flux density per beam of the source, set to *total_flux* (mJy/beam)
- *e_peak_flux*: beam corrected 1-sigma statistical error on the peak flux, set to *e_total_flux* (mJy/beam)
  of north, set to image *bpa* (degrees);
- *distance_from_center*: angular separation between the fitted position
  and the image pointing center (degrees)
- *polar_angle*: angle west of north of fitted position in image plane (degrees)
- *snr*: signal-to-noise ratio of the null detection; defined as
  *ave_total_flux* / *total_flux*



.. _error:

*****
error
*****
This is a look-up table containing explanations for each possible
*error_id* in the **image** table.

- *id*: referenced by the *error_id* column in the **image** table;
  updates cascade to that table
- *reason*: reason why an image was given that particular *error_id*

 id : reason                            

--------------------------------------------------------------------

-  1 : image missing necessary header keyword(s)
-  2 : number of visibilities < *min nvis*
-  3 : sensitivity metric (noise x sqrt(int. time)) <= 0 or > *max sensitivity metric*
-  4 : beam axis ratio > *max beam axis ratio*
-  5 : bad imaging target (NCP or planet)
-  6 : problem source in image field-of-view
-  7 : PyBDSF failed to process
-  8 : zero sources extracted
-  9 : source count metric > *max source count metric*
- 10 : number of CLEAN iterations < *min niter*
- 11 : image bmin in pixels < *min bpix* or > *max bpix* 
- 12 : image missing primary calibrators
- 13 : image missing CLEAN components
- 14 : image missing NX table


.. _image:

*****
image
*****
This table provides a record of every image processed by **vdp**.
The more useful keywords from the image header are summarized
in the table, as well.

- *id*: unique identifier for each new image; deletions cascade to
  the **detected_island** and **vlite_unique** tables
- *filename*: image filename with full directory path
- *imsize*: image size in pixels (pixels); header keywords
  (``NAXIS1``, ``NAXIS2``)
- *obs_ra*: right ascension of image pointing center (degrees);
  header keyword ``OBSRA``
- *obs_dec*: declination of image pointing center (degrees);
  header keyword ``OBSDEC``
- *glon*: galactic longitude of image pointing center (degrees);
  header keyword ``GLON`` or calculated if missing
- *glat*: galactic latitude of image pointing center (degrees);
  header keyword ``GLAT`` or calculated if missing
- *az_star*: azimuth of image pointing center at mjdtime (degrees);
  header keyword ``AZ_STAR``
- *el_star*: altitude of image pointing center at mjdtime (degrees);
  header keyword ``EL_STAR``
- *pa_star*: parallactic angle of image pointing center at mjdtime (degrees);
  header keyword ``PA_STAR``
- *az_end*: azimuth of image pointing center at end of observation (degrees);
  header keyword ``AZ_END``
- *el_end*: altitude of image pointing center at end of observation (degrees);
  header keyword ``EL_END``
- *pa_end*: parallactic angle of image pointing center at end of observation (degrees);
  header keyword ``PA_END``
- *az_i*: azimuth of image pointing center at mjdtime (degrees);
  calculated with astropy
- *alt_i*: altitude of image pointing center at mjdtime (degrees);
  calculated with astropy
- *parang_i*: parallactic angle of image pointing center at mjdtime (degrees);
  calculated with astropy
- *az_f*: azimuth of image pointing center at mjdtime+duration (degrees);
  calculated with astropy
- *alt_f*: altitude of image pointing center at mjdtime+duration (degrees);
  calculated with astropy
- *parang_f*: parallactic angle of image pointing center at mjdtime+duration (degrees);
  calculated with astropy
- *lst_i*: VLA local sidereal time at mjdtime (hrs);
  calculated with astropy
- *lst_f*: VLA local sidereal time at mjdtime+duration (hrs);
  calculated with astropy
- *pixel_scale*: number of arcseconds spanned by each pixel in the
  image (arcsec/pixel); header keyword ``CDELT1`` or ``CDELT2``
- *object*: name of the object being observed as given by the
  primary observer; header keyword ``OBJECT``
- *obs_date*: date observations were acquired formatted as
  yyyy-mm-dd; header keyword ``DATE-OBS``
- *map_date*: date image was created formatted as yyyy-mm-dd;
  header keyword ``DATE-MAP``
- *obs_freq*: rest frequency of the VLITE observations (MHz);
  header keyword ``RESTFREQ`` or ``CRVAL3`` or ``CRVAL4``
- *primary_freq*: frequency of the primary observations (GHz);
  taken from VLITE image filename
- *bmaj*: image beam FWHM major axis (arcsec); header keyword
  ``BMAJ`` or ``CLEANBMJ``
- *bmin*: image beam FWHM minor axis (arcsec); header keyword
  ``BMIN`` or ``CLEANBMN``
- *bpa*: image beam position angle measured east of north (degrees);
  header keyword ``BPA`` or ``CLEANBPA``
- *noise*: estimate of the rms noise measured in the center of the
  image (mJy/beam); header keyword ``ACTNOISE``
- *peak*: flux value of the brightest pixel in the image (mJy/beam);
  header keyword ``PEAK`` or ``DATAMAX``
- *config*: VLA configuration; header keyword ``CONFIG``
- *cycle*: VLITE configuration cycle, e.g. 1, 2, 3...
- *semester*: NRAO semester of the VLITE observation
- *nvis*: number of visibilities in the data before imaging; header
  keyword ``NVIS``
- *niter*: number of CLEAN iterations; header keyword ``NITER``
- *mjdtime*: Modified Julian Date at the start of the observations;
  Calculated from *obs_date* + header keyword ``STARTIME``
- *tau_time*: total integration time on source (sec); header
  keyword ``TAU_TIME``
- *duration*: total duration of the observations (sec); header
  keyword ``DURATION``
- *radius*: size of the circular field-of-view used in source finding
  (degrees); calculated as ((``NAXIS2`` / 2.) * *scale*) * ``CDELT2``,
  where *scale* is defined in the **pybdsf_params** section of the
  configuration file
- *nsrc*: number of sources found in the image by PyBDSF
- *nclean*: number of CLEANed sources in the image
- *rms_box*: size and step size of the box used by PyBDSF to estimate
  the image background mean and noise (pixels)
- *stage*: value of the highest **vdp** stage completed on the image:
  
  - 1 = reading the image
  - 2 = source finding
  - 3 = source association
  - 4 = catalog matching

- *catalogs_checked*: list of the names of other radio catalogs that
  have been checked for sources which can be positionally matched to
  the sources detected in this image
- *error_id*: value assigned if the image fails one of the quality
  checks; references the *id* column of the **error** table & is
  updated if that table is updated
- *nearest_problem*: name of the nearest source which is known
  to cause imaging problems for VLITE
- *separation*: angular separation between the VLITE image pointing
  center and the nearest known problem source (degrees)
- *pri_cals*: list of primary calibrators used
- *ass_flag*: True if image *bmin* within range allowed for source 
  association for image *config* and *cycle*
- *nsn*: Number of SN tables in UVOUT file. 
  Three SN tables mean A&P selfcal applied to image
- *tsky*: 341 MHz brightness temp [K] from the Global Sky Model 
  (de Oliveira-Costa et al. 2008)
- *sunsep*: angular separation between VLITE image pointing 
  center and the Sun at mjdtime (degrees)
- *pbkey*: Primary beam key for image
- *pb_flag*: True if primary beam calculation for image was possible
- *nvisnx*: Number of visibilities in the NX table; before self-cal applied
- *ninterval*: Number of observing intervals in NX table
- *max_dt*: Max time of observing intervals in NX table
- *nbeam*: Number of beams in primary beam image calculation
- *pbparangs*: Array of parallactic angles used in primary beam image calculation
- *pbweights*: Array of weighting factors used in primary beam image calculation

.. _run_config:

**********
run_config
**********
The contents of the configuration file are stored in this table
each time the pipeline is executed. Other information about the
run of the pipeline like the start time, execution time, and number
of images read are is also recorded.

- *id*: unique identifier corresponding to each run of **vdp**
- *config_file*: name of the configuration file used for the run
- *log_file*: name of the log file
- *start_time*: date and time the run was started; formatted as
  yyyy-mm-dd hh:mm:ss
- *execution_time*: length of the time taken to execute the full
  run; formatted as hh:mm:ss.s
- *nimages*: number of images read during the run
- *stages*: contents of the configuration file **stages** section
  stored as key-value pairs
- *options*: contents of the configuration file **options** section
  stored as key-value pairs
- *setup*: contents of the configuration file **setup** section
  stored as key-value pairs
- *pybdsf_params*: contents of the configuration file **pybdsf_params**
  section stored as key-value pairs
- *image_qa_params*:contents of the configuration file **image_qa_params**
  section stored as key-value pairs

.. _vlite_unique:

************
vlite_unique
************
Sources detected in VLITE images which are not successfully matched with
any sources in other radio catalogs are deemed to be "VLITE unique", or VU,
sources. These sources are placed into their own table to quickly isolate
transient candidates, steep spectrum, or other sources of potential interest.
A VLITE source is first added to this table when no matches are found during
the catalog matching stage. A new entry for the VU source is added to the
table for every image which could have or does contain the same VU source.

- *id*: unique row id number
- *image_id*: unique identifier for the VLITE image; references the *id*
  column of the **image** table & deletions cascade from that table
- *assoc_id*: references the *id* column of the **assoc_source** table
  in which the source has *nmatches* = 0; deletions cascade from that
  table
- *detected*: boolean ``True`` or ``False`` identifying whether or not
  the source was detected in the image
  
.. _sky_catalogs:

The "radcat" Schema
^^^^^^^^^^^^^^^^^^^
Cross-matching of VLITE sources with sources in other radio catalogs
and surveys is done within the database taking advantage of the Q3C
spatial indexing and functions. To enable this, all radio catalogs
used for comparison must be stored in the same database as the
VLITE sources. Since these catalogs are not part of the same
organizational structure as the main VLITE database, they are
stored in a separate schema with the name "radcat". Every
survey/catalog is a separate table with the same format in
this schema. There is also a table **radcat.catalogs** which
provides a list of every table contained in the schema along
with information about that catalog (telescope, frequency,
spatial resolution, publication reference).

.. _catalogs:

***************
radcat.catalogs
***************
This table contains information about every radio catalog
that is available for cross-matching with VLITE sources.
It can be used to get the name of the catalog referenced
by the *catalog_id* column in the **catalog_match** table.

- *id*: unique row id to identify the catalog
- *name*: name of the survey/catalog and table in the "radcat" schema
- *telescope*: name of the radio telescope used to acquire the
  catalog data
- *frequency*: frequency of the catalog observations (MHz)
- *resolution*: approximate spatial resolution of the catalog observations;
  used to ensure VLITE sources are only cross-matched with sources from
  radio catalogs in the same resolution class
- *reference*: reference to the publication containing the catalog

.. _catalog_sources:

*********************
radcat.[catalog_name]
*********************
Each survey/catalog included in the "radcat" schema is stored
in a separate table of the name **radcat.[catalog_name]**,
where **[catalog_name]** is the name of the catalog as it
appears in the **radcat.catalogs** table. Each one of these
tables has the same format with the following columns:

- *id*: row id number of the source referenced by the *src_id*
  column of the **catalog_match** table
- *name*: name given to the source
- *ra*: source right ascension (degrees)
- *e_ra*: error on the source right ascension (degrees)
- *dec*: source declination (degrees)
- *e_dec*: error on the source declination (degrees)
- *total_flux*: total integrated flux (mJy)
- *e_total_flux*: error on the total integrated flux (mJy)
- *peak_flux*: peak flux density per beam (mJy/beam)
- *e_peak_flux*: error on the peak flux density (mJy/beam)
- *maj*: size of the source semi-major axis (arcsec)
- *e_maj*: error on the semi-major axis size (arcsec)
- *min*: size of the source semi-minor axis (arcsec)
- *e_min*: error on the semi-minor axis size (arcsec)
- *pa*: source position angle (degrees)
- *e_pa*: error on the position angle
- *rms*: local noise estimate, or the rms noise in the image (mJy/beam)
- *field*: name of the field, or image, where the source was detected
- *catalog_id*: id number of the catalog that the source is in;
  references the *id* column of the **radcat.catalogs** table
- *pt_like*: True or false if catalog source is point-like or not. Null if unknown


