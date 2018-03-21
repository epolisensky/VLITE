# VLITE Database Pipeline (VDP)
### Scripts for automating source extraction and database storage.

#### Version 1.7

A post-processing pipeline to adaptively extract and catalog point sources is
being developed to enhance the scientific value and accessibility of data
products generated by the VLA Low-band Ionosphere and Transient Experiment
(VLITE) on the Karl G. Jansky Very Large Array (VLA). In contrast to other
radio sky surveys, the commensal observing mode of VLITE results in varying
depths, sensitivities, and spatial resolutions across the sky based on the
configuration of the VLA, location on the sky, and time on source specified
by the primary observer for their independent science objectives. Therefore,
previously developed tools and methods for generating source catalogs and
survey statistics are not always appropriate for VLITE's diverse and growing
set of data. A catalog of radio sources detected in every VLITE image
is created from source fit parameters measured using the Python Blob Detector
and Source Finder software (`PyBDSF`; Mohan & Rafferty 2015). Detected VLITE
sources are positionally associated with each other in a resolution-dependent
manner, and are cross-matched to other radio sky surveys to aid in the
detection of transient sources and enable creation of radio flux spectra for
sources across the entire northern sky.

![VLITE 3 years](https://github.com/erichards/VLITE/blob/develop/VLITE_3YEARS_map.png "VLITE 3 years")

## Files
This repository contains the following files under the vdp/ directory:

- **database/createdb.py**: Python functions to create the database tables
- **database/dbclasses.py**: Python module defining classes used in VDP and
their functions
- **database/dbio.py**: Python functions which handle database I/O for VDP
- **database/pybdsfcat.py**: Python function to read list of sources output
by `PyBDSF`
- **matching/matchfuncs.py**: Python module containing functions for positional
cross-matching
- **matching/radioxmatch.py**: Python driver script for performing association
and catalog cross-matching
- **skycatalog/catalogio.py**: Python module for reading other radio sky
survey source lists and catalogs
- **skycatalog/skycatdb.py**: Python module for creation of separate database
schema and insertion of other radio sky survey catalogs
- **sourcefinding/beam_tools.py**: Python functions to primary beam correct
VLITE flux measurements
- **sourcefinding/runbdsf.py**: Python script which runs `PyBDSF` source finding
- **OUTLINE.md**: outline of the VDP tasks
- **VDPDataFlow.jpg**: VDP data flow chart
- **VDPLogicFlow.jpg**: decision tree demonstrating VDP processing
stages & options
- **errors.py**: Python module containing class definitions for errors
- **example_config.yaml**: example `YAML` configuration file used to set
run parameters for VDP
- **maintenance.py**: Python module containing code to cluster & analyze
the indexed database tables.
- **vdp.py**: Python module that drives the whole pipeline
- **test_cmbranches.py**: Python testing script for catalog matching stage
- **test_sabranches.py**: Python testing script for source association stage
- **test_sfbranches.py**: Python testing script for source finding stage
- **timeout.py**: Python module defining decorator function to timeout long
running functions

## Running the Code
```
>> python vdp.py example_config.yaml
```
To see all optional command line arguments, type:
```
>> python vdp.py -h
```
or
```
>> python vdp.py --help
```

## Dependencies
- psycopg2
- yaml
- numpy
- astropy
- ephem
- pandas
- PyBDSF

