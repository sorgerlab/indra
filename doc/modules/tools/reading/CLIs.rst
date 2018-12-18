Command Line Interfaces for high-throughput reading
===================================================

Python CLI to run reading on local files
----------------------------------------

.. argparse::
    :module: indra.tools.reading.read_files
    :func: make_parser
    :prog: python -m indra.tools.reading.read_files

Python CLI to run the DRUM reading system
-----------------------------------------

.. argparse::
    :module: indra.tools.reading.run_drum_reading
    :func: make_parser
    :prog: python -m indra.tools.reading.run_drum_reading

Python CLI for submitting reading pipelines
-------------------------------------------

.. argparse::
    :module: indra.tools.reading.submit_reading_pipeline
    :func: create_parser
    :prog: python -m indra.tools.reading.submit_reading_pipeline

Python CLI to monitor running batch jobs
----------------------------------------

.. argparse::
    :module: indra.tools.reading.wait_for_complete
    :func: make_parser
    :prog: python -m indra.tools.reading.wait_for_complete


Python CLI to generate stats on reading results
-----------------------------------------------

.. argparse::
    :module: indra.tools.reading.util.reading_results_stats
    :func: make_parser
    :prog: python -m indra.tools.reading.util.reading_results_stats