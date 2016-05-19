import boto3
import zlib
from indra import reach
import pickle
import sys

pmid_list_file = sys.argv[1]

# Load the list of PMIDs from the given file
with open(pmid_list_file) as f:
    pmid_list = [line.strip('\n') for line in f.readlines()]

# Initialize S3 stuff
bucket_name ='bigmech'
client = boto3.client('s3')

stmts = {}
for pmid in pmid_list[0:10]:
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
    reach_proc = reach.process_json_str(reach_json, citation='PMID%s' % pmid)
    stmts[pmid] = reach_proc.statements

with open('reach_stmts.pkl', 'w') as f:
    pickle.dump(stmts, f)
