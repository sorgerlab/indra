import logging
import boto3

# Logger
logger = logging.getLogger('s3_client')

# Check Amazon credentials here?

# Create global boto3 client singleton? Perhaps by lazy initialization?

bucket_name ='bigmech'
client = boto3.client('s3')
prefix = 'papers/'

def check_pmid(pmid):
    if isinstance(pmid, int):
        pmid = str(pmid)
    if not pmid.startswith('PMID'):
        pmid = 'PMID' + str(pmid)
    return pmid

def get_full_text(pmid):
    pmid = check_pmid(pmid)
    fulltext_prefix = prefix + pmid + '/fulltext/'
    # Check for Open Access nxml source
    oa_xml_key = fulltext_prefix + 'pmc_oa_xml'
    oa_xml_s3obj = client.get_object(Bucket=bucket, Key=oa_xml_key)
    oa_xml_compressed = oa_xml_s3obj['Body'].read()
    # oa_xml = gzip_decompress(oa_xml_compressed)
    # Check if key not found. If not found, look for auth_xml
    return (oa_xml, 'xml')

def get_reach_output(pmid):
    pmid = check_pmid(pmid)
    reach_key = prefix + pmid + '/reach'
    reach_s3obj = client.get_object(Bucket=bucket, Key=reach_key)
    meta = reach_s3obj['Metadata']
    reach_result = reach_s3obj['Body'].read()
    return reach_result

def get_reach_version(pmid):
    pmid = check_pmid(pmid)
    reach_key = prefix + pmid + '/reach'
    try:
        reach_gz_obj = client.get_object(Key=reach_key, Bucket=bucket_name)
        logger.info("%s: found REACH output on S3; checking metadata" % pmid)
        reach_metadata = reach_gz_obj['Metadata']
        reach_version = reach_metadata.get('reach_version')
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            logger.info('No REACH output found on S3 for key' % reach_key)
            reach_version = None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    return reach_version


