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


Tools to run the DRUM reading system (:py:mod:`indra.tools.reading.run_drum_reading`)
-------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.run_drum_reading
    :members:

CLI to run the DRUM reading system
----------------------------------

.. argparse::
    :module: indra.tools.reading.run_drum_reading
    :func: make_parser
    :prog: python -m indra.tools.reading.run_drum_reading

Python tools for submitting reading pipelines (:py:mod:`indra.tools.reading.submit_reading_pipeline`)
-----------------------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.submit_reading_pipeline
    :members:

Python CLI for submitting reading pipelines
-------------------------------------------

.. argparse::
    :module: indra.tools.reading.submit_reading_pipeline
    :func: create_parser
    :prog: python -m indra.tools.reading.submit_reading_pipeline

Python cli to monitor running batch jobs
----------------------------------------

.. argparse::
    :module: indra.tools.reading.wait_for_complete
    :func: make_parser
    :prog: python -m indra.tools.reading.wait_for_complete
