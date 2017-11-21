# VLITE Post-Processing Pipeline (P3)
## For all our storage and analysis needs
### Task/Architecture Outline

*Note: Italics indicates functionality that has not yet been implemented.*

#### Stage 0: Pre-Run Checks
1. Read in [run configuration file](https://github.com/erichards/VLITE/blob/develop/p3/example_config.yaml)
2. Check if images exist/path is correct
3. Make sure requested cross-match catalogs exist
4. [Connect to database](https://docs.google.com/spreadsheets/d/e/2PACX-1vR20qGzJ7U3hFBNYZ1IUJWcFpdlOmfjQKv_8pk6aRW7BuljZ6VGNWyagHnsMVkZ6_Y9-Dl1vEwNv8Bg/pubhtml "database schema")
5. Create separate schema for sky catalogs and populate tables, if necessary

#### Stage 1: Image Cataloging, Organization, & *First Quality Check*
1. Initialize Image object
2. Check the Image table to see if the image has been processed before
3. If yes, check run configuration file to see if re-processing will be forced
4. If no, update database image table
5. *Perform quality checks and update database error table/image table error_id*

#### Stage 2: Source Finding & *Final Image Quality Check*
1. Initialize BDSFImage object with arguments for PyBDSF process_image task
read in from run configuration file
2. Run PyBDSF source finding
3. Convert PyBDSF class objects to DetectedSource objects
4. *Run quality checks based on the number density of detected sources*
5. Insert extracted sources into database detected_source and detected_island tables

#### Stage 3: Source Association
1. Set search area to same box used in source finding
2. Extract all previously detected VLITE sources from the database assoc_source
table that lie within the same region on the sky as the current image and come
from images with similar spatial resolutions
3. Cross-match current sources to the previously detected sources:
    1. If a match is found, update assoc_source table entry with *weighted*
    average for position and shape parameters, add one to the number of detections
    (ndetect), and update the detected_source table
    entry assoc_id to refer to the correct assoc_source row
    2. If no previous VLITE match is found, add the current source to the
    assoc_source table as a new entry and set ndetect to 1.

#### *Stage 4: Sky Survey Catalog Cross-Matching*
1. Cross-match every new detected source that was added to the assoc_source
table with every sky catalog table specified in the run configuration file
2. Extract sky catalog sources from the same sky region as the image (same
as Stage 3)
3. Add results to the database catalog_match table and update the number of
matches (nmatches) in the assoc_source table.