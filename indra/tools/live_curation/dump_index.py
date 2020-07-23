"""This is a script to dump all the corpora on S3 into an index file."""

import boto3

res = s3.list_objects(Bucket='world-modelers', Prefix='indra_models')
s3 = boto3.session.Session(profile_name='wm').client('s3')

corpora = []
for entry in res['Content']:
    if entry['Key'].endswith('/statements.json'):
        corpus_id = entry['Key'].split('/')[1]
        mod = entry['LastModified']
        corpora.append((corpus_id, mod.strftime('%Y-%m-%d-%H-%M-%S')))

with open('index.csv', 'w') as fh:
    for corpus_id, mod in sorted(corpora, key=lambda x: x[1]):
        fh.write('%s,%s\n' % (corpus_id, mod))
