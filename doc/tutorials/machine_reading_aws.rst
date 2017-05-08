Large-Scale Machine Reading with Amazon Batch
=============================================

The following doc describes the steps involved in reading a large numbers of
papers in parallel on Amazon EC2 using REACH, caching the JSON output on Amazon
S3, processing the REACH output into INDRA Statements, and then caching the
statements also on S3. Prerequisites for doing the following are:

* An Amazon S3 bucket containing full text contents for papers, keyed by
  Pubmed ID (creation of this S3 repository will be described in another
  tutorial).
* Amazon AWS credentials for using AWS Batch.
* A corpus of PMIDs (see *Large-Scale Machine Reading with Starcluster* for
  information on how to assemble this)
* Optional: Elsevier text and data mining API key and institution key for
  subscriber access to Elsevier full text content.

How it Works
------------

* The reading pipeline makes use of a Docker image that contains INDRA and all
  necessary dependencies, including REACH, Kappa, PySB, etc. The Docker file
  for this image is available at: https://github.com/johnbachman/indra_docker.
* The INDRA Docker image is built by AWS Codebuild and pushed to Amazon's
  EC2 Container Service (ECS), where it is available via the Repository URI::

    292075781285.dkr.ecr.us-east-1.amazonaws.com/indra

* An AWS Batch *Compute Environment* named "run_reach" is configured to use
  this Docker image for handling AWS jobs. This compute environment is configured
  to use only Spot instances with a maximum spot price of 40% of the on-demand
  price, and 16 vCPUs.

* An AWS *Job Queue*, "run_reach_queue", is configured to use instances of the
  "run_reach" Compute Environment.

* An AWS *Job Definition*, "run_reach_jobdef", is configured to run in the
  "run_reach_queue", and to use 16 vCPUs and 30GiB of RAM.

* Reading jobs are submitted by running the script:: 

    python -m indra.tools.reading.submit_reading_pipeline_aws read [args]

  which, given a list of PMIDs:

  * Copies the PMID list to the key ``reading_results/[job_name]/pmids`` on
    Amazon S3
  * Breaks the list up into chunks (e.g., of 3000 PMIDs) and submits an AWS
    Batch job for each (using the "run_reach_jobdef" definition as a template).

* The ECS instance created by the AWS Batch job runs the script
  ``indra.tools.reading.run_reach_on_pmids_aws``, which:

  * Checks for cached content on Amazon S3
  * If the PMID has not been read by the current version of REACH, checks for
    content
  * If the content is not available, downloads the content using the INDRA
    literature client, and caches on S3
  * The content to be read is downloaded to the ``/tmp`` directory of the
    instance
  * REACH is run using the command-line interface (RunReachCLI), and configured
    to read the papers in the ``/tmp`` directory using all of the vCPUs on the
    instance
  * When done, the result REACH JSON in the output folder is uploaded to S3
  * The JSON for both the previously and newly read papers is processed in
    parallel to INDRA Statements
  * The resulting subset of statements for the given range of papers is cached
    on S3 at ``reading_results/[job_name]/stmts/[start_ix]_[end_ix].pkl``. This
    set of statements takes the form of a pickled (protocol 3) Python dict
    with PMIDs as keys and lists of INDRA Statements as values.
  * In addition, information about the sources of content available for each
    PMID is cached for each PMID subset at
    ``reading_results/[job_name]/content_types/[start_ix]_[end_ix].pkl``.

* When the reading jobs for each of the subsets of PMIDs have been completed
  and cached on S3, the final combined set of statements (and combined
  information on content sources) can be assembled using::

    python -m indra.tools.reading.submit_reading_pipeline_aws combine [job_name]

  * This script submits an AWS batch job for a machine with 1 vCPU but a large
    amount of memory (60GiB)
  * The job runs the script ``indra.tools.reading.assemble_reach_stmts_aws``,
    which unpickles the results from all of the PMID subsets, combines them,
    and stores them on S3
  * The resulting files are obtainable from S3 at
    ``reading_results/[job_name]/stmts.pkl`` and
    ``reading_results/[job_name]/content_types.pkl``.

* To run the entire pipeline, where the assembly of the combined set of
  statements is automatically performed after the reading step is completed,
  run::

    python -m indra.tools.reading.submit_reading_pipeline_aws full [args]

