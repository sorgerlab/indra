Large-Scale Machine Reading
===========================

The following doc describes the steps involved in reading a large numbers of
papers in parallel on Amazon EC2 using REACH, caching the JSON output on Amazon
S3, then processing the REACH output into INDRA Statements. Prerequisites for
doing the following are:

* A cluster of Amazon EC2 nodes configured using Starcluster, with INDRA
  installed and in the PYTHONPATH
* An Amazon S3 bucket containing full text contents for papers, keyed by
  Pubmed ID (creation of this S3 repository will be described in another
  tutorial).

This tutorial goes through the individual steps involved before describing how
all of them can be run through the use of a single submission script,
submit_reading_pipeline.py.

Note also that the prerequisite installation steps can be streamlined by
putting them in a setup script that can be re-run upon instantiating a new
Amazon cluster or by using them to configure a custom Amazon EC2 AMI.

Install REACH
-------------

Install SBT. On an EC2 Linux machine, run the following lines (drawn from
http://www.scala-sbt.org/0.13/docs/Installing-sbt-on-Linux.html)::

    echo "deb https://dl.bintray.com/sbt/debian /" | sudo tee -a /etc/apt/sources.list.d/sbt.list
    sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 642AC823
    sudo apt-get update
    sudo apt-get install sbt

Clone REACH from https://github.com/clulab/reach. Under reach/project, create a
file, assembly.sbt, containing the following line::

    addSbtPlugin("com.eed3si9n" % "sbt-assembly" % "0.14.3")

Add the following lines to reach/build.sbt::

    test in assembly := {}
    mainClass in assembly := Some("org.clulab.reach.ReachCLI")

This prevents Scala from running all the tests during the SBT assembly step,
and assigns ReachCLI as the main class.

Compile and assemble REACH. Note that the path to the .ivy2 directory must be
given. Use the assembly task to assemble a fat JAR containing all of the
dependencies with the correct main class. Run the following from the directory
containing the REACH build.sbt file (e.g., /pmc/reach).::

    sbt -Dsbt.ivy.home=/pmc/reach/.ivy2 compile
    sbt -Dsbt.ivy.home=/pmc/reach/.ivy2 assembly

Install Amazon S3 support
-------------------------

Install boto3::

    pip install boto3

.. note::

    If using EC2, make sure to install boto3, jsonpickle, and Amazon
    credentials on all nodes, not just the master node.

Add Amazon credentials to access the S3 bucket. First create the .aws directory
on the EC2 instance::

    mkdir /home/sgeadmin/.aws

Then set up Amazon credentials, for example by copying from your local machine
using StarCluster::

    starcluster put mycluster ~/.aws/credentials /home/sgeadmin/.aws

Install other dependencies
--------------------------

::

    pip install jsonpickle # Necessary to process JSON from S3
    pip install --upgrade jnius-indra # Necessary for REACH
    export JAVA_HOME=/usr/lib/jvm/java-1.7.0-openjdk-amd64

Assemble a Corpus of PMIDs
--------------------------

The first step in large-scale reading is to put together a file containing
relevant Pubmed IDs. The simplest way to do this is to use the Pubmed search
API to find papers associated with particular gene names, biological processes,
or other search terms.

For example, to assemble a list of papers for SOS2 curated in Entrez Gene
that are available in the Pubmed Central Open Access subset:

.. ipython:: python

    from indra.literature import *

    # Pick an example gene
    gene = 'SOS2'

    # Get a list of PMIDs for the gene
    pmids = pubmed_client.get_ids_for_gene(gene)

    # Get the PMIDs that have XML in PMC
    pmids_oa_xml = pmc_client.filter_pmids(pmids, 'oa_xml')

    # Write the results to a file
    with open('%s_pmids.txt' % gene, 'w') as f:
        for pmid in pmids_oa_xml:
            f.write('%s\n' % pmid)

This creates a file, SOS2_pmids.txt, containing the PMIDs that we will read
with REACH.

Process the papers with REACH
-----------------------------

The next step is to read the content of the papers with REACH in a
parallelizable, high-throughput way. To do this, run the script
indra/tools/reading/run_reach_on_pmids.py. If necessary update the lines at the
top of the script with the REACH settings, e.g.::

    cleanup = False
    verbose = True
    path_to_reach = '/pmc/reach/target/scala-2.11/reach-assembly-1.3.2-SNAPSHOT.jar'
    reach_version = '1.3.2'
    source_text = 'pmc_oa_xml'
    force_read = False

The reach_version is important because it is used to determine whether the
paper has already been read with this version of REACH (in which case it will
be skipped), or if the REACH output needs to be updated. Alternatively, if you
want to read all the papers regardless of whether they've been read before with
the given version of REACH, set the force_read variable to True.

Next, create a top-level temporary directory to use during reading. This will
be used to store the input files and the JSON output::

    mkdir my_temp_dir

Run run_reach_on_pmids.py, passing arguments for the PMID list file, the temp
directory, the number of cores to use on the machine, the PMID start index (in
the PMID list file) and the end index. The start and end indices are used to
subdivide the job into parallelizable chunks. If the end index is greater than
the total number of PMIDs, it will process up to the last one in the list. For
example::

    python run_reach_on_pmids.py SOS2_pmids.txt my_temp_dir 8 0 10

This uses 8 cores to process the first ten papers listed in the file
SOS2_pmids.txt. REACH will run, output the JSON files in the temporary
directory, e.g. in my_temp_dir/read_0_to_10_MSP6YI/output, assemble the JSON
files together, and upload the results to S3. If you attempt to process the
files again with the same version of REACH, the script will detect that the
JSON output from that version is already on S3 and skip those papers.

This can be submitted to run offline using the job scheduler on EC2 with, e.g.::

    qsub -b y -cwd -V -pe orte 8 python run_reach_on_pmids.py SOS2_pmids.txt my_temp_dir 8 0 10

.. note::

    The number of cores requested in the qsub call ('-pe orte 8') should match
    the number of cores passed to the run_reach_on_pmids.py script, which
    determines the number of threads that REACH will attempt to use (the
    third-to-last argument above). This should also match the total number of
    nodes on the Amazon EC2 node (e.g., 8 cores for c3.2xlarge). This way the
    job scheduler will schedule the job to run on all the cores of a single EC2
    node, and REACH will use them all.

Extract INDRA Statements from the REACH output on S3
----------------------------------------------------

The script indra/tools/reading/process_reach_from_s3.py is used to extract
INDRA Statements from the REACH output uploaded to S3 in the previous step.
This process can also be parallelized by submitting chunks of papers to be
processed by different cores. The INDRA statements for each chunk of papers are
pickled and can be assembled into a single pickle file in a subsequent step.

Following the example above, run the following to process the REACH output
for the SOS2 papers into INDRA statements. We'll do this in two chunks to
show how the process can be parallelized and the statements assembled from
multiple files::

    python process_reach_from_s3.py SOS2_pmids.txt 0 5
    python process_reach_from_s3.py SOS2_pmids.txt 5 10

The two runs create two different files for the results from the seven papers,
reach_stmts_0_5.pkl (with statements from the first five papers) and
reach_stmts_5_7.pkl (with statements from the last two). Note that the results
are pickled as a dict (rather than a list), with PMIDs as keys and lists of
Statements as values.

Of course, what we really want is a single file containing all of the
statements for the entire corpus. To get this, run::

    python assemble_reach_stmts.py reach_stmts_*.pkl

The results will be stored in reach_stmts.pkl.

Running the whole pipeline with one script
------------------------------------------

If you want to run the whole pipeline in one go, you can run the script
submit_reading_pipeline.py (in indra/tools/reading) on a cluster of Amazon
EC2 nodes. The script divides up the jobs evenly among the nodes and cores.
Usage::

    python submit_reading_pipeline.py pmid_list tmp_dir num_nodes num_cores_per_node

For example if you have a cluster with 8 c3.8xlarge nodes with 32 VCPUs each,
you would call it with::

    python submit_reading_pipeline.py SOS2_pmids.txt my_tmp_dir 8 32

The script submits the jobs to the scheduler with appropriate dependencies
such that the REACH reading step completes first, then the INDRA processing
step, and then the final assembly into a single pickle file.


