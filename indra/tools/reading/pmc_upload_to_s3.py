from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

if __name__ == '__main__':
    import gzip
    import boto3
    import sys
    import pymysql
    import time
    from io import StringIO

    if len(sys.argv) != 3:
        print("Usage: %s start_index end_index_inclusive" % sys.argv[0])
        sys.exit(1)

    client = boto3.client('s3')
    db = pymysql.connect(host='master', user='root', db='pmc')
    cursor = db.cursor(pymysql.cursors.DictCursor)

    bucket_name = 'bigmech'

    start_index = int(sys.argv[1])
    end_index = int(sys.argv[2])

    cursor.execute("select * from papers where id >= %d and id <= %d and pmid != ''" %
                    (start_index, end_index))

    start = time.time()
    for rec in cursor.fetchall():
        my_id = rec['id']
        if int(my_id) % 10 == 0:
            print('%s %s' % (my_id, rec['pmid']))
        # Shortcut for making key names
        def key_prefix(s):
            return 'papers/PMID%s/%s' % (rec['pmid'], s)
        # Store PMCID
        pmcid_key = key_prefix('pmcid')
        client.put_object(Key=pmcid_key, Body=rec['pmcid'], Bucket=bucket_name)
        # Store DOI
        doi = rec['doi']
        if doi:
            doi_key = key_prefix('doi')
            client.put_object(Key=doi_key, Body=doi, Bucket=bucket_name)
        # Store full text 
        with open(rec['path']) as f:
            xml = f.read()
        buf = StringIO()
        content_type = rec['content_type']
        gzf = gzip.GzipFile(content_type, 'wb', 6, buf)
        gzf.write(xml)
        gzf.close()
        xml_key = key_prefix('fulltext/%s' % content_type)
        client.put_object(Key=xml_key, Body=buf.getvalue(), Bucket=bucket_name)
    end = time.time()
    elapsed = end - start
    print("Elapsed: %s" % elapsed)

