Welcome to the VLITE Database Pipeline Documentation
====================================================
The VLITE Database Pipeline, or **vdp**, is a collection of
Python scripts to automate the measurement and database
archiving of radio-emitting astronomical point sources detected in
VLITE images.

**vdp** uses the Python Blob Detector and Source Finder
(`PyBDSF <http://www.astron.nl/citt/pybdsf/>`_;
Mohan & Rafferty 2015) for source finding and PostgreSQL
for database storage. The sky indexing scheme `Q3C
<https://github.com/segasai/q3c>`_ (Koposov, S., &
Bartunov, O. 2006) implemented for PostgreSQL is
utilized to enable high performance queries and
positional matching across tables as the VLITE catalog
grows.

This documentation provides a high level overview of how
to run **vdp**, as well as detailed descriptions
of each stage in the pipeline.

For more information about VLITE in general, see `here.
<http://vlite.nrao.edu/>`_

Quick Start
-----------

Execution of the pipeline is controlled through a YAML
configuration file. Here's an example file to get you
started: :download:`example_config.yaml <example_config.yaml>`

The pipeline is started from the command line with the
configuration file::

  $ python vdp.py example_config.yaml

Check out the `VDPGuide IPython notebook
<https://github.com/erichards/VLITE/blob/master/vdp/notebooks/VDPGuide.ipynb>`_
for a step-by-step walkthrough of the pipeline stages and
analysis examples.


Getting Started with the VLITE Database Pipeline
================================================

.. toctree::
   :maxdepth: 3
	      
   basic_usage
   stages
   database
   usage_examples
   sql
   instructions

   
Modules
=======

.. toctree::
   :maxdepth: 2

   modules


Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

