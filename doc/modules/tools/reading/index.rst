High-throughput reading tools (:py:mod:`indra.tools.reading`)
==============================================================

INDRA defines interfaces to many text reading tools, however many of those only
handle reading at small scales. These tools are developed to harness reading at
arbitrary scales.

Tools used to run reading on a set of locally stored files (:py:mod:`indra.tools.reading.read_files`)
-----------------------------------------------------------------------------------------------------

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


Python tools for submitting reading pipelines (:py:mod:`indra.tools.reading.submit_reading_pipeline`)
-----------------------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.submit_reading_pipeline
    :members:



A pipeline that uses AWS Batch and caches on S3 (:py:mod:`indra.tools.reading.pmid_reading`)
--------------------------------------------------------------------------------------------

This pipeline makes use of AWS Batch jobs to scale readings arbitrarily, and
optimizes the reading by caching results on S3, thus preventing the user from
needing to reread content unnecessarily. This pipeline may someday be retired
in favor of the RDS reading pipeline (see below), however at this time this
method is nominally maintained.

The machinery for reading a list of PMIDs with REACH or SPARSER (:py:mod:`indra.tools.reading.pmid_reading.read_pmids`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: indra.tools.reading.pmid_reading.read_pmids
    :members:


A pipeline that uses AWS Batch and RDS (:py:mod:`indra_db.reading`)
-------------------------------------------------------------------

This pipeline no longer exists under the umbrella of INDRA, but rather lives in
the INDRA DB repo:

    https://github.com/indralab/indra_db/tree/master/indra_db/reading

More information can be found in the `indra_db` documentation. Unlike the other
pipelines, this system is aimed at continuous automatic reading of all content.


A pipeline that uses StarCluster (:py:mod:`indra.tools.reading.starcluster_reading`)
------------------------------------------------------------------------------------

The first pipeline developed for large scale reading. This method is no longer
actively maintained.


Low-level utilities used for or by reading pipelines
====================================================

General tools to help with bulk reading.


Tools to analyze logs (:py:mod:`indra.tools.reading.util.log_analysis_tools`)
-----------------------------------------------------------------------------

.. automodule:: indra.tools.reading.util.log_analysis_tools
    :members:


Tools to export statements by reach rule (:py:mod:`indra.tools.reading.util.export_stmts_by_reach_rule`)
--------------------------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.util.export_stmts_by_reach_rule
    :members:


Class for generating reports on a Batch reading job (:py:mod:`indra.tools.reading.util.reporter`)
-------------------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.util.reporter
    :members:


Tools used in creating CLI scripts for reading (:py:mod:`indra.tools.reading.util.script_tools`)
------------------------------------------------------------------------------------------------

.. automodule:: indra.tools.reading.util.script_tools
    :members:
