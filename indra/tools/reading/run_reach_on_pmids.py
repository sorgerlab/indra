# Iterate over the PMIDs
# Check in S3
# If present, copy to NXML folder
# Prepare temporary .conf file using appropriate number of cores
# Run reach on NXML files, send to output folder
# Join the resulting .json files
# Upload the .json to S3, mark down in a folder as being run

import os
import sys
import zlib
import tempfile
import shutil
import boto3
import botocore
import subprocess
import glob
import json
import cStringIO
import gzip

cleanup = False
verbose = True
path_to_reach = '/pmc/reach/target/scala-2.11/reach-assembly-1.3.2-SNAPSHOT.jar'
reach_version = '1.3.2'
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
# Initialize S3 stuff
bucket_name ='bigmech'
client = boto3.client('s3')
pmids_to_read = []

# Check if we've read the PMIDs already
for pmid in pmids_in_range:
    # See if we've already read this one
    reach_key = 'papers/PMID%s/reach' % pmid
    try:
        reach_gz_obj = client.get_object(Key=reach_key, Bucket=bucket_name)
        print "Found the object, so let's check the metadata"
        reach_metadata = reach_gz_obj['Metadata']
        if reach_metadata.get('reach_version') is not None and \
           reach_metadata.get('reach_version') == reach_version:
            print "Same version as what we've got, so don't do anything"
            continue
        # No version info, or not the same as the current version
        else:
            print "Not the same version, so read it!"
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            print "No object found for key %s, read it!" % reach_key
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    pmids_to_read.append(pmid)

if not pmids_to_read:
    print "No pmids to read!"
    sys.exit(0)

# Now iterate over the pmids to read  and download from S3 to the input
# directory
for pmid in pmids_to_read:
    # Look for the full text
    key_prefix = 'papers/PMID%s/fulltext' % pmid
    key_name = '%s/%s' % (key_prefix, source_text)
    print key_name
    # Get the content from S3
    try:
        xml_gz_obj = client.get_object(Key=key_name, Bucket=bucket_name)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            print "No object found for key %s" % key_name
            continue
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    # Get the content from the object
    xml_gz = xml_gz_obj['Body'].read()
    # Decode the gzipped content
    xml = zlib.decompress(xml_gz, 16+zlib.MAX_WBITS)
    # Write the contents to a file
    xml_path = os.path.join(input_dir, 'PMID%s.nxml' % pmid)
    with open(xml_path, 'w') as f:
        f.write(xml)

# Create the REACH configuration file
conf_file_text = """
# this is the directory that stores the raw nxml files
# this directory *must* exist
nxmlDir = {input_dir}
# this is where the output files containing the extracted mentions will be
# stored
# if this directory doesn't exist it will be created
friesDir = {output_dir}
# this is where the context files will be stored
# if this directory doesn't exist it will be created
contextDir = {output_dir}
# this is where the brat standoff and text files are dumped
bratDir = {output_dir}
# the encoding of input and output files
encoding = "utf-8"
# nxml2fries configuration
nxml2fries {{
  # this is a list of sections that we should ignore
  ignoreSections = ["references", "materials", "materials|methods", "methods", "supplementary-material"]
}}
# context engine config
contextEngine {{
    type = Policy4
    params = {{
        bound = 3
    }}
}}
# the output format for mentions: text, fries, indexcard (default is 'text')
outputType = "fries"
# this log file gets overwritten every time ReachCLI is executed
# so you should copy it if you want to keep it around
logFile = {base_dir}/log.txt
# grounding configuration
grounding: {{
  # List of AdHoc grounding files to insert, in order, into the grounding
  # search sequence. Each element of the list is a map of KB filename and
  # optional meta info (not yet used):
  #   example: {{ kb: "adhoc.tsv", source: "NMZ at CMU" }}
  adHocFiles: [
#     {{ kb: "adhoc.tsv", source: "NMZ at CMU" }}
  ]
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
        print '@@', line
(p_out, p_err) = p.communicate()
if p.returncode:
    raise Exception(p_out + '\n' + p_err)

# At this point, we have a directory full of JSON files
# Collect all the prefixes into a set, then iterate over the prefixes

def join_parts(prefix):
    """Join different REACH output JSON files into a single JSON."""
    entities = json.load(open(prefix + '.uaz.entities.json'))
    events = json.load(open(prefix + '.uaz.events.json'))
    sentences = json.load(open(prefix + '.uaz.sentences.json'))
    full = {'events': events, 'entities': entities, 'sentences': sentences}
    return full

def gzip_string(content, name):
    buf = cStringIO.StringIO()
    gzf = gzip.GzipFile(name, 'wb', 6, buf)
    gzf.write(content)
    gzf.close()
    return buf.getvalue()

# Collect prefixes
json_files = glob.glob(os.path.join(output_dir, '*.json'))
json_prefixes = set([])
for json_file in json_files:
    filename = os.path.basename(json_file)
    prefix = filename.split('.')[0]
    json_prefixes.add(prefix)
# Now iterate over the collected prefixes, combine the JSON, and send to S3
for json_prefix in json_prefixes:
    prefix_with_path = os.path.join(output_dir, json_prefix)
    full_json = join_parts(prefix_with_path)
    full_json_gz = gzip_string(json.dumps(full_json), 'reach_output.json')
    reach_key = 'papers/%s/reach' % json_prefix
    reach_metadata = {'reach_version': reach_version,
                      'source_text': source_text}
    client.put_object(Key=reach_key, Body=full_json_gz, Bucket=bucket_name,
                      Metadata=reach_metadata)

if cleanup:
    shutil.rmtree(base_dir)


