import logging
import boto3
import zlib
import botocore
import json
import cStringIO
import gzip

# Logger
logger = logging.getLogger('s3_client')

# Check Amazon credentials here?

# Create global boto3 client singleton? Perhaps by lazy initialization?

bucket_name ='bigmech'
client = boto3.client('s3')
prefix = 'papers/'
s3 = boto3.resource('s3')
bucket = s3.Bucket(bucket_name)

def check_pmid(pmid):
    if isinstance(pmid, int):
        pmid = str(pmid)
    if not pmid.startswith('PMID'):
        pmid = 'PMID' + str(pmid)
    return pmid


def get_pmid_key(pmid):
    pmid = check_pmid(pmid)
    return prefix + pmid


def filter_keys(prefix):
    return list(bucket.objects.filter(Prefix=prefix))


def check_key(key):
    try:
        s3.Object(bucket_name, key).load()
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == '404':
            exists = False
        else:
            raise e
    else:
        exists = True
    return exists


def get_full_text(pmid, full_text_type='pmc_oa_xml'):
    pmid = check_pmid(pmid)
    oa_xml_key = prefix + pmid + '/fulltext/' + full_text_type
    # Check for Open Access nxml source
    try:
        xml_gz_obj = client.get_object(Bucket=bucket_name, Key=oa_xml_key)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            logger.info('get_full_text: no object found for key %s' %
                        oa_xml_key)
            return None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    # Get the content from the object
    xml_gz = xml_gz_obj['Body'].read()
    # Decode the gzipped content
    xml = zlib.decompress(xml_gz, 16+zlib.MAX_WBITS)
    return xml.decode('utf8')


def put_full_text(pmid, text, full_text_type='pmc_oa_xml'):
    pmid = check_pmid(pmid)
    xml_key = prefix + pmid + '/fulltext/' + full_text_type
    xml_gz = gzip_string(text, '%s.nxml' % pmid)
    client.put_object(Key=xml_key, Body=xml_gz, Bucket=bucket_name)


def put_abstract(pmid, text):
    xml_key = get_pmid_key(pmid) + '/abstract'
    xml_gz = gzip_string(text, '%s.nxml' % pmid)
    client.put_object(Key=xml_key, Body=xml_gz, Bucket=bucket_name)


def get_reach_version(pmid):
    pmid = check_pmid(pmid)
    reach_key = prefix + pmid + '/reach'
    try:
        reach_gz_obj = client.get_object(Key=reach_key, Bucket=bucket_name)
        logger.info("%s: found REACH output on S3; checking version" % pmid)
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


def get_reach_output(pmid):
    pmid = check_pmid(pmid)
    reach_key = prefix + pmid + '/reach'
    reach_s3obj = client.get_object(Bucket=bucket, Key=reach_key)
    meta = reach_s3obj['Metadata']
    reach_result = reach_s3obj['Body'].read()
    return reach_result


def put_reach_output(reach_output, pmid_key, reach_version, source_text):
    full_json_gz = gzip_string(json.dumps(reach_output), 'reach_output.json')
    reach_key = 'papers/%s/reach' % pmid_key
    reach_metadata = {'reach_version': reach_version,
                      'source_text': source_text}
    client.put_object(Key=reach_key, Body=full_json_gz, Bucket=bucket_name,
                      Metadata=reach_metadata)


def gzip_string(content, name):
    buf = cStringIO.StringIO()
    gzf = gzip.GzipFile(name, 'wb', 6, buf)
    gzf.write(content.encode('utf8'))
    gzf.close()
    return buf.getvalue()


