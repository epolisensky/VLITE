# VLITE
Scripts for post-processing and analyzing VLITE images.

## Post-Processing Pipeline (P3)
See p3driver.py.

### Stage 1: Source Finding
Source finding executed using Python Blob Detector and Source Finder (PyBDSF; github.com/lofar-astron/PyBDSF) -- SourceFinding/runPyBDSF.py

### Stage 2: Database Insertion
Source finding results are stored in an SQLite database -- SourceFinding/database.py

### Stage 3: Image Quality Check
Filter images based on quality -- TBD

### Stage 4: Catalog Cross-Matching
Cross-match detected sources to other radio frequency sky survey catalogs -- CatalogMatching/radioXmatch.py

