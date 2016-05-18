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

cleanup = True

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
pmids_to_read = pmid_list[start_index:end_index]
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
# Iterate over the pmids and download from S3 to the input directory
for pmid in pmids_to_read:
    key_prefix = 'papers/PMID%s/fulltext' % pmid
    content_type = 'pmc_oa_xml'
    key_name = '%s/%s' % (key_prefix, content_type)
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

#if cleanup:
#    shutil.rmtree(base_dir)

