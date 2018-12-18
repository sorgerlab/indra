High-throughput reading CLI's
=============================

The reading tools include several Python CLI's to run reading or other related
tasks.


Run reading on local files
--------------------------

.. argparse::
    :module: indra.tools.reading.read_files
    :func: make_parser
    :prog: python -m indra.tools.reading.read_files

Run REACH and/or SPARSER locally on a list of PMIDs using S3 caching
--------------------------------------------------------------------

.. argparse::
    :module: indra.tools.reading.pmid_reading.read_pmids
    :func: make_parser
    :prog: python -m indra.tools.reading.pmid_reading.read_pmids


Submit AWS Batch reading jobs
-----------------------------

.. argparse::
    :module: indra.tools.reading.submit_reading_pipeline
    :func: create_parser
    :prog: python -m indra.tools.reading.submit_reading_pipeline


Monitor running AWS Batch jobs
------------------------------

.. argparse::
    :module: indra.tools.reading.wait_for_complete
    :func: make_parser
    :prog: python -m indra.tools.reading.wait_for_complete


Run the DRUM reading system
---------------------------

.. argparse::
    :module: indra.tools.reading.run_drum_reading
    :func: make_parser
    :prog: python -m indra.tools.reading.run_drum_reading


Generate stats on AWS Batch reading results
-------------------------------------------

.. argparse::
    :module: indra.tools.reading.util.reading_results_stats
    :func: make_parser
    :prog: python -m indra.tools.reading.util.reading_results_stats
