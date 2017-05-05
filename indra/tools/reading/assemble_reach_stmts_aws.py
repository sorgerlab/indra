
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

if __name__ == '__main__':
    from indra.tools.reading import run_reach_on_pmids as rr
    import boto3
    import botocore
    import os
    import sys
    import pickle
    import logging

    client = boto3.client('s3')
    bucket_name = 'bigmech'
    basename = sys.argv[1]

    # The trailing slash here is important
    prefix = 'reading_results/%s/stmts/' % basename
    # Get all keys associated with reading results
    stmt_file_keys = []
    marker = ''
    while True:
        res = client.list_objects(Bucket=bucket_name, Prefix=prefix,
                                  Delimiter='/', Marker=marker)
        if not res.get('Contents'):
            break
        keys = [item['Key'] for item in res['Contents']]
        stmt_file_keys += keys
        # If the response is truncated, get another round of results
        if not res.get('IsTruncated'):
            marker = res['NextMarker']
        else:
            break
    # Now that we have the keys, get and unpickle each of the pickles
    stmts = {}
    for key in stmt_file_keys:
        try:
            stmt_obj = client.get_object(Bucket=bucket_name, Key=key)
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] =='NoSuchKey':
                logger.debug('key %s not in S3' % key)
                continue
            # If there was some other kind of problem, re-raise the exception
            else:
                raise e
        stmt_bytes = stmt_obj['Body'].read()
        stmt_subset = pickle.loads(stmt_bytes)
        # Add stmts into master list
        stmts.update(stmt_subset)

    # Write out the final statement set
    # Pickle the statements to a bytestring
    if stmts:
        pickle_key_name = 'reading_results/%s/stmts.pkl' % basename
        stmts_bytes = pickle.dumps(stmts)
        client.put_object(Key=pickle_key_name, Body=stmts_bytes,
                          Bucket=bucket_name)
