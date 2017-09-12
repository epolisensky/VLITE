# VLITE Post-Processing Pipeline (P3)
## For all our storage and analysis needs
### Task/Architecture Outline

*Note*: Italics indicates functionality that has not yet been implemented.

##### Pre-requisites:
- *run configuration file*
- SkyCatalogs database -- separate table for each sky survey catalog
- VLITE database -- all observations at all configurations

#### Stage 0: Pre-Run Checks
1. *Read in run configuration file*
2. Check path to database/create it if doesn't exist yet
3. Check if images exist/path is correct
4. Check if path to SkyCatalogs database is correct
5. Make sure requested cross-match catalogs exist in the SkyCatalogs database

#### Stage 1: Image Cataloging & Organization
1. Initialize Image object
2. Check the Image table to see if the image has been processed before
3. *If yes, check run configuration file to see if re-processing will be forced*
4. If no, update database Image table with image parameters

#### *Stage 2: Preliminary Image Quality Check*
1. Check nvis, # of antennas, tau_time, etc. to remove images with too little
integration time
2. Check for proximity to bright radio sources
3. Noise check (how?)
4. Update database Error table/Image table error_id
5. Possible future development: use this stage to figure out optimal PyBDSF
input parameters

#### Stage 3: Source Finding
1. Initialize BDSFImage object with arguments for PyBDSF process_image task
*read in from run configuration file*
2. *Choose processing mode specified in run config. file (i.e. try multiple
rms_box sizes, process residual image after original, etc.)*
3. Call PyBDSF
4. Convert PyBDSF class objects to DetectedSource objects
5. Insert detected sources into database rawSource and rawIsland tables

#### *Stage 4: Final Image Quality Check*
1. Run some tests based on the number density of detected sources? -- Just
looking for extreme outliers here that may have been missed by the first
quality checks. We can remove this step if experience demonstrates that this
is useless or too difficult to define generally.
2. Other ideas?
3. Update database Error table/Image table error_id
4. Possible future development: use this stage to figure out optimal matching
parameters. 

#### *Stage 5: Source Association*
1. Calculate radius of field-of-view for SQL cone search.
2. Extract all previously detected VLITE sources within the search range
and which were detected in the same configuration as the current sources from
the AssocSource table.
3. Cross-match current sources to extracted AssocSource sources:
    1. If a match is found, compute weighted average for position and shape,
    add (src_id, image_id) to list of rawSource_ids to link back to rawSource
    table for flux measurements, and add one to num_detections.
    2. If no previous VLITE match is found, add the current source to the
    AssocSource table as a new entry and set num_detections to 1.
4. Add one to num_nulls column for all un-matched AssocSource sources. This is
a source that was detected previously in an image covering the same area on the
sky in the same configuration, but was not detected in the current image.
Flagging these will help identify false detections and transients, assuming we
figure out how to properly account for different sensitivities.
5. Future development will hopefully include a way to cross-match sources
detected in different configurations, using the highest res. measurements as
the "true" source position and shape. Some thought and experience is needed to
determine how to best handle the resulting one-to-many (i.e. lower-to-higher
res.), many-to-one (i.e. higher-to-lower res.), and many-to-many associations
that will result when trying to do this. There are current implementation
examples for dealing with this, but they aren't great.

#### Stage 6: Sky Survey Catalog Cross-Matching
1. *Consult run configuration file to determine if catalog cross-matching is to
be run for all detected sources or just new sources.*
2. Extract catalog sources from the first sky survey table specified in the run
configuration file using same range in RA & Dec spanned by image used to
extract sources from the AssocSource table in Stage 5.
2. Cross-match VLITE sources with the extracted catalog sources.
3. Update the catalog_id, match_id, min_DeRuiter columns in AssocSource table
with match results.


### Pipeline Output
1. Updated database with the following tables:
    * Image -- Image metadata and record of reason why image was not able to be
    processed, if applicable.
    * rawIsland -- Flux, background, and noise measurements of islands of
    emission in an image from which sources are fit in PyBDSF. Multiple sources
    may come from a single island in a single image. The primary key is the
    combination of isl_id and image_id, where image_id is a foreign key linking
    to the Image table.
    * rawSource -- Fitting results from PyBDSF source finding. Each row
    corresponds to a single source detected in a single image, so this table
    may contain "duplicates" of the same source detected in multiple images.
    The isl_id is a foreign key linking to the rawIsland table. The primary key
    is the combination of src_id and image_id, where image_id is a foreign key
    linking to the Image table.
    * AssocSource -- Amalgamation of source measurements where multiple
    detections of the same source (in the same configuration) are combined to
    form a single position and shape. Each detected source is cross-matched
    with sources in this table to determine whether the source has been
    detected previously in the same configuration. The number of detections for
    a source (in each configuration) is recorded. The number of null detections
    for a source is also recorded (see Stage 5). Each source has a list of
    row_id foreign keys which link back to the corresponding rawSource
    measurements for extraction of lighcurves. Results of the sky survey
    catalog cross-matching are included as the foreign keys catalog_id and
    match_id, which can be used link to the corresponding matched catalog
    source in the SkyCatalogs database, and the minimum de Ruiter radius
    calculated between the associated sources.
    * Error -- List of reasons why an image could be rejected for processing.
    The error_id foreign key in the Image table links to the primary key id in
    this table.
2. PyBDSF outputs:
    * ASCII text file containing the list of sources from each image
    * ds9 regions file of the detected sources
    * residual image after subtracting the fitted sources
    * log file for each image
3. ds9 regions file of the matched catalog sources

*Additional optional outputs:*
1. ds9 regions file of the un-matched detected sources
2. model image of fitted sources from PyBDSF