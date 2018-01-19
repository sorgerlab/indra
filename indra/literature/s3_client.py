from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import logging
import boto3
import zlib
import botocore
import json
import gzip
from io import BytesIO
from xml.etree import ElementTree as ET
from indra import literature as lit
from indra.literature import elsevier_client, pubmed_client
from indra.util import UnicodeXMLTreeBuilder as UTB
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
    # As of boto3 1.4.4, there appears to be an incompatibility with the
    # future newstr string type (which is meant to be unicode-compatible
    # in Python 2) and the urllib quote function, which boto3 uses to
    # URL encode keys when making requests to S3. See, for example this
    # related error:
    # http://stackoverflow.com/questions/15115588/urllib-quote-throws-keyerror
    # In Python 2, this can be solved by explicitly convert to either to
    # unicode or to a Python 2 bytestring (bytes) explicitly.
    if sys.version_info.major < 3:
        pmid = unicode(pmid)
    return pmid


def get_pmid_key(pmid):
    pmid = check_pmid(pmid)
    return prefix + pmid


def get_reader_key(reader, pmid):
    return get_pmid_key(pmid) + '/' + reader


def filter_keys(prefix):
    return list(bucket.objects.filter(Prefix=prefix))


def get_upload_content(pmid, force_fulltext_lookup=False):
    """Get full text and/or abstract for paper and upload to S3."""
    # Make sure that the PMID doesn't start with PMID so that it doesn't
    # screw up the literature clients
    if pmid.startswith('PMID'):
        pmid = pmid[4:]
    # First, check S3:
    (ft_content_s3, ft_content_type_s3) = get_full_text(pmid)
    # The abstract is on S3 but there is no full text; if we're not forcing
    # fulltext lookup, then we're done
    if ft_content_type_s3 == 'abstract' and not force_fulltext_lookup:
        return (ft_content_s3, ft_content_type_s3)
    # If there's nothing (even an abstract on S3), or if there's an abstract
    # and we're forcing fulltext lookup, do the lookup
    elif ft_content_type_s3 is None or \
            (ft_content_type_s3 == 'abstract' and force_fulltext_lookup) or \
            (ft_content_type_s3 == 'elsevier_xml' and
                    not elsevier_client.extract_text(ft_content_s3)):
        if ft_content_type_s3 == 'elsevier_xml':
            logger.info('PMID%s: elsevier_xml cached on S3 is missing full '
                        'text element, getting again.' % pmid)
        # Try to retrieve from literature client
        logger.info("PMID%s: getting content using literature client" % pmid)
        (ft_content, ft_content_type) = lit.get_full_text(pmid, 'pmid')
        assert ft_content_type in ('pmc_oa_xml', 'elsevier_xml',
                                   'abstract', None)
        # If we tried to get the full text and didn't even get the abstract,
        # then there was probably a problem with the web service. Try to
        # get the abstract instead:
        if ft_content_type is None:
            return (None, None)
        # If we got the abstract, and we already had the abstract on S3, then
        # do nothing
        elif ft_content_type == 'abstract' and ft_content_type_s3 == 'abstract':
            logger.info("PMID%s: found abstract but already had it on " \
                        "S3; skipping" % pmid)
            return (ft_content, ft_content_type)
        # If we got the abstract, and we had nothing on S3, then upload
        elif ft_content_type == 'abstract' and ft_content_type_s3 is None:
            logger.info("PMID%s: found abstract, uploading to S3" % pmid)
            put_abstract(pmid, ft_content)
            return (ft_content, ft_content_type)
        # If we got elsevier_xml, but cannot get a full text element, then
        # get and put the abstract
        elif ft_content_type == 'elsevier_xml' and \
                not elsevier_client.extract_text(ft_content):
            logger.info("PMID%s: Couldn't get a full text element for "
                        "the elsevier_xml content; getting abstract "
                        % pmid)
            abstract = pubmed_client.get_abstract(pmid)
            # Abstract is None, so return None
            if abstract is None:
                logger.info("PMID%s: Unable to get abstract, returning None"
                            % pmid)
                return (None, None)
            # Otherwise, upload and return the abstract
            else:
                logger.info("PMID%s: Uploading and returning abstract "
                            % pmid)
                put_abstract(pmid, abstract)
                return (abstract, 'abstract')
        # We got a viable full text
        # (or something other than None or abstract...)
        else:
            logger.info("PMID%s: uploading and returning %s"
                        % (pmid, ft_content_type))
            put_full_text(pmid, ft_content, full_text_type=ft_content_type)
            return (ft_content, ft_content_type)
    # Some form of full text is already on S3
    else:
        # TODO
        # In future, could check for abstract even if full text is found, and
        # upload it just to have it
        return (ft_content_s3, ft_content_type_s3)
    # We should always return before we get here
    assert False


def get_gz_object(key):
    try:
        gz_obj = client.get_object(Bucket=bucket_name, Key=key)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            logger.debug('key %s not in S3' % key)
            return None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    # Get the content from the object
    gz_body = gz_obj['Body'].read()
    # Decode the gzipped content
    content = zlib.decompress(gz_body, 16+zlib.MAX_WBITS)
    return content.decode('utf8')


def get_object_metadata(key):
    try:
        metadata = client.head_object(Bucket=bucket_name, Key=key)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] =='NoSuchKey':
            logger.debug('key %s not in S3' % key)
            return None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    return metadata


def get_full_text(pmid, metadata=False):
    pmid = check_pmid(pmid)
    # Check for Open Access nxml source
    ft_prefix = get_pmid_key(pmid) + '/fulltext/'
    ft_objs = filter_keys(ft_prefix)
    # We have at least one full text
    if len(ft_objs) > 0:
        ft_keys = [ft_obj.key for ft_obj in ft_objs]
        # Look for full texts in order of desirability
        for content_type in ('pmc_oa_xml', 'pmc_auth_xml', 'elsevier_xml',
                             'pmc_oa_txt', 'txt'):
            ft_key = ft_prefix + content_type
            # We don't have this type of full text, move on
            if ft_key not in ft_keys:
                continue
            # We have this type of full text, so get it and return
            else:
                if not metadata:
                    content = get_gz_object(ft_key)
                else:
                    content = get_object_metadata(ft_key)
                if content:
                    logger.info('%s: found %s on S3' % (pmid, content_type))
                    return (content, content_type)
                else:
                    logger.info('%s: error getting %s' %
                                (pmid, content_type))
                    return (None, None)
        # If we've gotten here, it means there were full text keys not
        # included in the above
        logger.error('Unrecognized full text key %s for %s' %
                    (ft_keys, pmid))
        return (None, None)
    else:
        logger.debug('%s: no full texts found on S3, trying abstract' % pmid)
        abstract_key = get_pmid_key(pmid) + '/abstract'
        if not metadata:
            abstract = get_gz_object(abstract_key)
        else:
            abstract = get_object_metadata(abstract_key)
        if abstract is None:
            logger.info('%s: no full text or abstract found on S3' % pmid)
            return (None, None)
        else:
            logger.info('%s: found abstract on S3' % pmid)
            return (abstract, 'abstract')


def put_full_text(pmid, text, full_text_type='pmc_oa_xml'):
    pmid = check_pmid(pmid)
    xml_key = prefix + pmid + '/fulltext/' + full_text_type
    xml_gz = gzip_string(text, '%s.nxml' % pmid) # Encodes to UTF-8
    client.put_object(Key=xml_key, Body=xml_gz, Bucket=bucket_name)


def put_abstract(pmid, text):
    xml_key = get_pmid_key(pmid) + '/abstract'
    xml_gz = gzip_string(text, '%s.nxml' % pmid)
    client.put_object(Key=xml_key, Body=xml_gz, Bucket=bucket_name)


def get_reader_metadata(reader, pmid):
    reading_key = get_reader_key(reader, pmid)
    try:
        reading_gz_obj = client.get_object(Key=reading_key, Bucket=bucket_name)
        reading_metadata = reading_gz_obj['Metadata']
        # The reader version string comes back as str in Python 2, not unicode
        # Using str (instead of .decode) should work in both Python 2 and 3
        reader_version = reading_metadata.get('reader_version')
        logger.info("%s: found %s output on S3 for version %s" %
                    (pmid, reader, reader_version))
        if reader_version is not None:
            reader_version = str(reader_version)
        source_text = reading_metadata.get('source_text')
        if source_text is not None:
            source_text = str(source_text)
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == 'NoSuchKey':
            logger.info('No %s output found on S3 for key %s'
                        % (reader, reading_key))
            reader_version = None
            source_text = None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    return (reader_version, source_text)


def get_reader_output(reader, pmid):
    # Get the REACH JSON as unicode
    reader_json_str = get_reader_json_str(reader, pmid)
    # Now create the JSON object--the resulting obj will contain un-escaped
    # unicode data
    if reader_json_str is None:
        return None
    else:
        reader_json = json.loads(reader_json_str)
        return reader_json


def get_reader_json_str(reader, pmid):
    reader_key = get_reader_key(reader, pmid)
    try:
        reader_s3obj = client.get_object(Bucket=bucket_name, Key=reader_key)
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == 'NoSuchKey':
            logger.info('No %s output found on S3 for key %s'
                        % (reader, reader_key))
            return None
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    #meta = reader_s3obj['Metadata']
    reader_gz = reader_s3obj['Body'].read()
    # Gunzip the the content
    reader_bytes = zlib.decompress(reader_gz, 16+zlib.MAX_WBITS)
    # Convert from bytes to str (shouldn't affect content since all
    # Unicode should be escaped in the JSON)
    reader_uni = reader_bytes.decode('utf-8')
    return reader_uni


def put_reader_output(reader, reader_output, pmid, reader_version, source_text):
    if not isinstance(reader_version, basestring):
        raise ValueError("Reader version must be a string.")
    full_json_gz = gzip_string(json.dumps(reader_output),
                               '%s_output.json' % reader)
    reading_key = get_reader_key(reader, pmid)
    reading_metadata = {'reader_version': reader_version,
                        'source_text': source_text}
    client.put_object(Key=reading_key, Body=full_json_gz, Bucket=bucket_name,
                      Metadata=reading_metadata)


def gzip_string(content, name):
    buf = BytesIO()
    gzf = gzip.GzipFile(name, 'wb', 6, buf)
    gzf.write(content.encode('utf8'))
    gzf.close()
    return buf.getvalue()


