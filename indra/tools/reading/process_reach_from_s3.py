import boto3
import botocore
import zlib
from indra import reach
import pickle
import sys

usage = "Usage: %s pmid_list start_index end_index" % sys.argv[0]
if len(sys.argv) < 4:
    print usage
    sys.exit(1)

pmid_list_file = sys.argv[1]
start_ix = int(sys.argv[2])
end_ix = int(sys.argv[3])

# Load the list of PMIDs from the given file
with open(pmid_list_file) as f:
    pmid_list = [line.strip('\n') for line in f.readlines()]
if end_ix > len(pmid_list):
    end_ix = len(pmid_list)

# Initialize S3 stuff
bucket_name ='bigmech'
client = boto3.client('s3')

stmts = {}
for ix, pmid in enumerate(pmid_list[start_ix:end_ix]):
    # Get the reach output
    reach_key = 'papers/PMID%s/reach' % pmid
    try:
        print "Downloading reach output for %s" % pmid
        reach_gz_obj = client.get_object(Key=reach_key, Bucket=bucket_name)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            print "Reach output not found for PMID %s, skipping" % pmid
            continue
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    # Get the content from the object
    reach_gz = reach_gz_obj['Body'].read()
    # Decode the gzipped content
    reach_json = zlib.decompress(reach_gz, 16+zlib.MAX_WBITS)
    print "Processing"
    try:
        reach_proc = reach.process_json_str(reach_json,
                                            citation='PMID%s' % pmid)
    # If there's a problem, skip it
    except Exception as e:
        print "Exception processing %s" % pmid
        print e
        continue

    stmts[pmid] = reach_proc.statements

with open('reach_stmts_%d_%d.pkl' % (start_ix, end_ix), 'w') as f:
    pickle.dump(stmts, f)
