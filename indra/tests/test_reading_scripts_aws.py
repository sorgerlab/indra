
import boto3
from os import path
from subprocess import check_call
from indra.db import util as dbu
from indra.util import zip_string


s3 = boto3.client('s3')


HERE = path.dirname(path.abspath(__file__))


def test_db_reading_help():
    check_call(['python', '-m', 'indra.tools.reading.db_reading.read_db_aws',
                '--help'])


def test_db_reading_noraml_querty():
    # Put some basic stuff in the test databsae
    N = 6
    db = dbu.get_test_db()
    db._clear(force=True)
    db.copy('text_ref', [(i, 'PMID80945%d' % i) for i in range(N)],
            cols=('id', 'pmid'))
    text_content = [
        (i, i, 'pubmed', 'text', 'abstract',
         zip_string('MEK phosphorylates ERK in test %d.' % i))
        for i in range(N)
        ]
    db.copy('text_content', text_content,
            cols=('id', 'text_ref_id', 'source', 'format', 'text_type',
                  'content'))

    # Put an id file on s3
    basename = 'local_test_run'
    s3_inp_prefix = 'reading_inputs/%s/' % basename
    s3_out_prefix = 'reading_results/%s/' % basename
    local_dir = path.join(HERE, 'local_test')
    s3.put_object(Bucket='bigmech', Key=s3_inp_prefix + 'id_list',
                  Body='\n'.join(['tcid: %d' % i for i in range(5)]))

    # Call the reading tool
    check_call(['python', '-m', 'indra.tools.reading.db_reading.read_db_aws',
                basename, local_dir, 'unread', 'all', '4', '0', '2', '-r',
                'sparser', '--test'])

    # Remove garbage on s3
    for pref in [s3_inp_prefix, s3_out_prefix]:
        res = s3.list_objects(Bucket='bigmech', Prefix=pref)
        for entry in res['Contents']:
            print("Removing %s..." % entry['Key'])
            s3.delete_object(Bucket='bigmech', Key=entry['Key'])
