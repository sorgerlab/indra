High-throughput reading tools (:py:mod:`indra.tools.reading`)
==============================================================

High-throughput reading pipelines
---------------------------------

.. toctree::
   :maxdepth: 2

   pmid_reading/index
   starcluster_reading/index


Reading utilities
-----------------

.. toctree::
   :maxdepth: 2

   util/index

Run reading on a set of locally stored files (:py:mod:`indra.tools.reading.read_files`)
---------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.read_files
    :members:

Classes defining and implementing interfaces to different readers (:py:mod:`indra.tools.reading.readers`)
---------------------------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.readers
    :members:


Run the DRUM reading system (:py:mod:`indra.tools.reading.run_drum_reading`)
----------------------------------------------------------------------------

.. automodule:: indra.tools.reading.run_drum_reading
    :members:

Python tools for submitting reading pipelines (:py:mod:`indra.tools.reading.submit_reading_pipeline`)
-----------------------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.submit_reading_pipeline
    :members:

Python cli for submitting reading pipelines (:py:mod:`indra.tools.reading.submit_reading_pipeline`)
---------------------------------------------------------------------------------------------------

.. argparse::
    :module: indra.tools.reading.submit_reading_pipeline
    :func: create_parser
    :prog: python submit_reading_pipeline.py

Python cli to monitor running batch jobs (:py:mod:`indra.tools.reading.wait_for_complete`)
------------------------------------------------------------------------------------------

.. argparse::
    :module: indra.tools.reading.wait_for_complete
    :func: make_parser
    :prog: python wait_for_complete.py
