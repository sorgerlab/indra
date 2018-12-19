Remote Reading Pipelines
========================

There are three reading pipelines that have been developed for reading on
remote high-performance systems.


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


