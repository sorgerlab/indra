# Iterate over the PMIDs
# Check in S3
# If present, copy to NXML folder
# Prepare temporary .conf file using appropriate number of cores
# Run reach on NXML files, send to output folder
# Join the resulting .json files
# Upload the .json to S3, mark down in a folder as being run

import os
import sys
import tempfile
import shutil
import boto3
import botocore
import subprocess
import glob
import json
import logging
from indra.literature import pmc_client, s3_client, get_full_text

cleanup = False
verbose = True
#path_to_reach = '/pmc/reach/target/scala-2.11/reach-assembly-1.3.2-SNAPSHOT.jar'
path_to_reach = '/Users/johnbachman/Dropbox/1johndata/Knowledge File/Biology/Research/Big Mechanism/reach/target/scala-2.11/reach-assembly-1.3.2-SNAPSHOT.jar'
reach_version = '1.3.2'
force_read = True
source_text = 'pmc_oa_xml'

# Check the arguments
usage = "Usage: %s pmid_list tmp_dir num_cores start_index end_index" \
         % sys.argv[0]
if len(sys.argv) < 6:
    print usage
    sys.exit(1)
# Get the command line arguments
(pmid_list_file, tmp_dir, num_cores, start_index, end_index) = sys.argv[1:]
start_index = int(start_index)
end_index = int(end_index)
num_cores = int(num_cores)

# Logger
logger = logging.getLogger('runreach')

# Load the list of PMIDs from the given file
with open(pmid_list_file) as f:
    pmid_list = [line.strip('\n') for line in f.readlines()]
if end_index > len(pmid_list):
    end_index = len(pmid_list)
pmids_in_range = pmid_list[start_index:end_index]

# Create the temp directories for input and output
base_dir = tempfile.mkdtemp(prefix='read_%s_to_%s_' % (start_index, end_index),
                            dir=tmp_dir)
input_dir = os.path.join(base_dir, 'input')
output_dir = os.path.join(base_dir, 'output')
os.makedirs(input_dir)
os.makedirs(output_dir)

pmids_to_read = []

# If we're re-reading no matter what, we don't have to check for existing
# REACH output
if force_read:
    pmids_to_read = pmids_in_range
# Otherwise, check if we've read the PMIDs already
else:
    for pmid in pmids_in_range:
        found_reach_version = s3_client.get_reach_version(pmid)
        # Found it, same version
        if found_reach_version is not None and \
           found_reach_version == reach_version:
            logger.info('%s: found same version (%s), skipping' %
                        (pmid, found_reach_version))
        # Found it, different version
        else:
            logger.info('%s: found %s, current %s; will re-read' %
                        (pmid, found_reach_version, reach_version))
            pmids_to_read.append(pmid)

if not pmids_to_read:
    logger.info('No pmids to read!')
    sys.exit(0)
else:
    logger.info('Preparing to run REACH on %d PMIDs' % len(pmids_to_read))

# Now iterate over the pmids to read and download from S3 to the input
# directory
num_found_s3 = 0
num_found_not_s3 = 0
for pmid in pmids_to_read:
    # Look for the full text
    xml = s3_client.get_full_text(pmid)
    # If we don't find the XML on S3, look for it using the PMC client
    if xml:
        num_found_s3 += 1
    else:
        logger.info('No content for %s from S3' % pmid)
        (content, content_type) = get_full_text(pmid, 'pmid')
        if content_type == 'nxml':
            logger.info('Found nxml for %s from PMC web service' % pmid)
            xml = content
            num_found_not_s3 += 1
            # Upload the xml to S3 for next time
            logger.info('Uploading full text for %s to S3' % pmid)
            s3_client.put_full_text(pmid, xml, full_text_type='pmc_oa_xml')
        else:
            logger.info('No full text found for %s' % pmid)
    # Write the contents to a file
    if xml:
        xml_path = os.path.join(input_dir, 'PMID%s.nxml' % pmid)
        with open(xml_path, 'w') as f:
            # The XML string is Unicode
            enc = xml.encode('utf8')
            f.write(enc)

logger.info('Found full text for %d PMIDs (%d S3, %d other)' %
            ((num_found_s3 + num_found_not_s3), num_found_s3,
              num_found_not_s3))

# Create the REACH configuration file
conf_file_text = """
#
# Configuration file for reach
#

# this is the directory that stores the raw nxml, .csv, and/or .tsv files
# this directory *must* exist
papersDir = {input_dir}

# this is where the output files containing the extracted mentions will be stored
# if this directory doesn't exist it will be created
outDir = {output_dir}

# the output format for mentions: text, fries, indexcard, or assembly-csv (default is 'fries')
outputType = "fries"

# whether or not assembly should be run
withAssembly = false

# this is where the context files will be stored
# if this directory doesn't exist it will be created
contextDir = {output_dir}

# this is where the brat standoff and text files are dumped
bratDir = {output_dir}

# verbose logging
verbose = true

# the encoding of input and output files
encoding = "utf-8"


# this is a list of sections that we should ignore
ignoreSections = ["references", "materials", "materials|methods", "methods", "supplementary-material"]

# context engine config
contextEngine {{
    type = Policy4
    params = {{
        bound = 3
    }}
}}

# this log file gets overwritten every time ReachCLI is executed
# so you should copy it if you want to keep it around
logFile = {base_dir}/log.txt

# grounding configuration
grounding: {{
  # List of AdHoc grounding files to insert, in order, into the grounding search sequence.
  # Each element of the list is a map of KB filename and optional meta info (not yet used):
  #   example: {{ kb: "adhoc.tsv", source: "NMZ at CMU" }}
  adHocFiles: [
#    {{ kb: "NER-Grounding-Override.tsv.gz", source: "MITRE/NMZ/BG feedback overrides" }}
  ]

  # flag to turn off the influence of species on grounding
  overrideSpecies = true
}}

# number of simultaneous threads to use for parallelization
threadLimit = {num_cores}

# ReadPapers
ReadPapers.papersDir = src/test/resources/inputs/nxml/
ReadPapers.serializedPapers = mentions.ser
""".format(base_dir=base_dir, input_dir=input_dir, output_dir=output_dir,
           num_cores=num_cores)

# Write the configuration file to the temp directory
conf_file_path = os.path.join(base_dir, 'indra.conf')
with open(conf_file_path, 'w') as f:
    f.write(conf_file_text)

# Run REACH!
args = ['java', '-jar', path_to_reach, conf_file_path]
p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if verbose:
    for line in iter(p.stdout.readline, b''):
        logger.info(line)
(p_out, p_err) = p.communicate()
if p.returncode:
    raise Exception(p_out + '\n' + p_err)

# At this point, we have a directory full of JSON files
# Collect all the prefixes into a set, then iterate over the prefixes
def join_parts(prefix):
    """Join different REACH output JSON files into a single JSON."""
    try:
        entities = json.load(open(prefix + '.uaz.entities.json'))
        events = json.load(open(prefix + '.uaz.events.json'))
        sentences = json.load(open(prefix + '.uaz.sentences.json'))
    except IOError as e:
        logging.error('Failed to open JSON files for %s; REACH error?' % prefix)
        return None
    return {'events': events, 'entities': entities, 'sentences': sentences}

# Collect prefixes
json_files = glob.glob(os.path.join(output_dir, '*.json'))
json_prefixes = set([])
for json_file in json_files:
    filename = os.path.basename(json_file)
    prefix = filename.split('.')[0]
    json_prefixes.add(prefix)
# Now iterate over the collected prefixes, combine the JSON, and send to S3
num_uploaded = 0
for json_prefix in json_prefixes:
    prefix_with_path = os.path.join(output_dir, json_prefix)
    full_json = join_parts(prefix_with_path)
    if full_json is None:
        num_failures += 1
    else:
        s3_client.put_reach_output(full_json, json_prefix, reach_version,
                                   source_text)
        num_uploaded += 1

logger.info('Uploaded REACH JSON for %d files to S3 (%d failures)' %
            (num_uploaded, num_failures))

if cleanup:
    shutil.rmtree(base_dir)


