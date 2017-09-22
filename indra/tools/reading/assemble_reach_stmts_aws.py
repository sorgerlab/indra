
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str


def assemble_batch_results(result_type):
    # The trailing slash here is important
    prefix = 'reading_results/%s/%s/' % (basename, result_type)
    # Get all keys associated with reading results
    result_file_keys = []
    marker = ''
    # Page through the files
    while True:
        res = client.list_objects(Bucket=bucket_name, Prefix=prefix,
                                  Delimiter='/', Marker=marker)
        if not res.get('Contents'):
            break
        result_file_keys += [item['Key'] for item in res['Contents']]
        # If the response is truncated, get another round of results
        if res.get('IsTruncated'):
            marker = res['NextMarker']
        else:
            break
    # Now that we have the keys, get and unpickle each of the pickles
    results = {}
    for key in result_file_keys:
        try:
            logger.info('Downloading and unpickling %s' % key)
            result_obj = client.get_object(Bucket=bucket_name, Key=key)
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] =='NoSuchKey':
                logger.debug('Key %s not found on S3' % key)
                continue
            # If there was some other kind of problem, re-raise the exception
            else:
                raise e
        result_bytes = result_obj['Body'].read()
        result_subset = pickle.loads(result_bytes)
        # Add results into master list
        results.update(result_subset)

    # Write out the final statement set
    # Pickle the statements to a bytestring
    if results:
        pickle_key_name = 'reading_results/%s/%s.pkl' % (basename, result_type)
        logger.info('Pickling combined file %s' % pickle_key_name)
        results_bytes = pickle.dumps(results)
        logger.info('Uploading combined file %s' % pickle_key_name)
        client.put_object(Key=pickle_key_name, Body=results_bytes,
                          Bucket=bucket_name)


if __name__ == '__main__':
    import boto3
    import botocore
    import sys
    import pickle
    import logging

    logger = logging.getLogger('assemble_reach')

    client = boto3.client('s3')
    bucket_name = 'bigmech'
    basename = sys.argv[1]

    result_types = ('content_types', 'stmts')
    for rt in result_types:
        assemble_batch_results(rt)
