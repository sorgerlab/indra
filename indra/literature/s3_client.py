from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import boto3
import zlib
import botocore
import json
import gzip
from io import BytesIO

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

# Logger
logger = logging.getLogger('s3_client')

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


def get_reach_key(pmid):
    return get_pmid_key(pmid) + '/reach'


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


def get_gz_object(key):
    try:
        gz_obj = client.get_object(Bucket=bucket_name, Key=key)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            logger.info('get_full_text: no object found for key %s' %
                        key)
            return None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    # Get the content from the object
    gz_body = gz_obj['Body'].read()
    # Decode the gzipped content
    content = zlib.decompress(gz_body, 16+zlib.MAX_WBITS)
    return content.decode('utf8')


def get_full_text(pmid):
    pmid = check_pmid(pmid)
    # Check for Open Access nxml source
    ft_prefix = get_pmid_key(pmid) + '/fulltext/'
    ft_objs = filter_keys(ft_prefix)
    # We have at least one full text
    if len(ft_objs) > 0:
        ft_keys = [ft_obj.key for ft_obj in ft_objs]
        # Look for full texts in order of desirability
        for content_type in ('pmc_oa_xml', 'pmc_auth_xml', 'pmc_oa_txt',
                             'txt'):
            ft_key = ft_prefix + content_type
            # We don't have this type of full text, move on
            if ft_key not in ft_keys:
                continue
            # We have this type of full text, so get it and return
            else:
                content = get_gz_object(ft_key)
                if content:
                    logger.info('Found %s for %s' % (content_type, pmid))
                    return (content, content_type)
                else:
                    logger.info('Error getting %s for %s' %
                                (content_type, pmid))
                    return (None, None)
        # If we've gotten here, it means there were full text keys not
        # included in the above
        logger.info('Unrecognized full text key %s for %s' %
                    (ft_keys, pmid))
        return (None, None)
    else:
        logger.info('No full texts found for %s, trying abstract' % pmid)
        abstract_key = get_pmid_key(pmid) + '/abstract'
        abstract = get_gz_object(abstract_key)
        if abstract is None:
            logger.info('No abstract found for %s' % pmid)
            return (None, None)
        else:
            logger.info('Found abstract for %s' % pmid)
            return (abstract, 'abstract')


def put_full_text(pmid, text, full_text_type='pmc_oa_xml'):
    pmid = check_pmid(pmid)
    xml_key = prefix + pmid + '/fulltext/' + full_text_type
    xml_gz = gzip_string(text, '%s.nxml' % pmid)
    client.put_object(Key=xml_key, Body=xml_gz, Bucket=bucket_name)


def put_abstract(pmid, text):
    xml_key = get_pmid_key(pmid) + '/abstract'
    xml_gz = gzip_string(text, '%s.nxml' % pmid)
    client.put_object(Key=xml_key, Body=xml_gz, Bucket=bucket_name)


def get_reach_metadata(pmid):
    reach_key = get_reach_key(pmid)
    try:
        reach_gz_obj = client.get_object(Key=reach_key, Bucket=bucket_name)
        logger.info("%s: found REACH output on S3; checking version" % pmid)
        reach_metadata = reach_gz_obj['Metadata']
        # The REACH version string comes back as str in Python 2, not unicode
        # Using str (instead of .decode) should work in both Python 2 and 3
        reach_version = reach_metadata.get('reach_version')
        if reach_version is not None:
            reach_version = str(reach_version)
        source_text = reach_metadata.get('source_text')
        if source_text is not None:
            source_text = str(source_text)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            logger.info('No REACH output found on S3 for key' % reach_key)
            reach_version = None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    return (reach_version, source_text)


def get_reach_output(pmid):
    # Get the REACH JSON as unicode
    reach_json_str = get_reach_json_str(pmid)
    # Now create the JSON object--the resulting obj will contain un-escaped
    # unicode data
    if reach_json_str is None:
        return None
    else:
        reach_json = json.loads(reach_json_str)
        return reach_json


def get_reach_json_str(pmid):
    reach_key = get_reach_key(pmid)
    try:
        reach_s3obj = client.get_object(Bucket=bucket_name, Key=reach_key)
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            logger.info('No REACH output found on S3 for key' % reach_key)
            return None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    #meta = reach_s3obj['Metadata']
    reach_gz = reach_s3obj['Body'].read()
    # Gunzip the the content
    reach_bytes = zlib.decompress(reach_gz, 16+zlib.MAX_WBITS)
    # Convert from bytes to str (shouldn't affect content since all
    # Unicode should be escaped in the JSON)
    reach_uni = reach_bytes.decode('utf-8')
    return reach_uni


def put_reach_output(reach_output, pmid, reach_version, source_text):
    if not isinstance(reach_version, basestring):
        raise ValueError("REACH version must be a string.")
    full_json_gz = gzip_string(json.dumps(reach_output), 'reach_output.json')
    reach_key = get_reach_key(pmid)
    reach_metadata = {'reach_version': reach_version,
                      'source_text': source_text}
    client.put_object(Key=reach_key, Body=full_json_gz, Bucket=bucket_name,
                      Metadata=reach_metadata)


def gzip_string(content, name):
    buf = BytesIO()
    gzf = gzip.GzipFile(name, 'wb', 6, buf)
    gzf.write(content.encode('utf8'))
    gzf.close()
    return buf.getvalue()


