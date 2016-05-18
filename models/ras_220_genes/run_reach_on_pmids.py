# Iterate over the PMIDs
# Check in S3
# If present, copy to NXML folder
# Prepare temporary .conf file using appropriate number of cores
# Run reach on NXML files, send to output folder
# Join the resulting .json files
# Upload the .json to S3, mark down in a folder as being run

import sys
import boto3
import botocore

# Check the arguments
usage = "Usage: %s pmid_list txt_dir output_dir num_cores start_index end_index" \
         % sys.argv[0]
if len(sys.argv) < 7:
    print usage
    sys.exit(1)
# Get the command line arguments
(pmid_list_file, txt_dir, output_dir,
                 num_cores, start_index, end_index) = sys.argv[1:]
start_index = int(start_index)
end_index = int(end_index)
num_cores = int(num_cores)
# Load the list of PMIDs from the given file
with open(pmid_list_file) as f:
    pmid_list = [line.strip('\n') for line in f.readlines()]
pmids_to_read = pmid_list[start_index:end_index]
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
    xml_gz = xml_gz['Body'].read()
    # Decode the gzipped content
    buf = cStringIO.StringIO()

    # Have we already read this paper? With this version of REACH?
    # If not, do nothing

