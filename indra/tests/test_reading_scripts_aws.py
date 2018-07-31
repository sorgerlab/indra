
import boto3
from os import path, chdir
from subprocess import check_call
from nose.plugins.attrib import attr

from indra.db import util as dbu
from indra.util import zip_string
from indra.tools.reading import submit_reading_pipeline as srp


s3 = boto3.client('s3')


HERE = path.dirname(path.abspath(__file__))


@attr('nonpublic')
def test_db_reading_help():
    chdir(path.expanduser('~'))
    check_call(['python', '-m', 'indra.tools.reading.db_reading.read_db_aws',
                '--help'])


@attr('nonpublic')
def test_normal_db_reading_call():
    chdir(path.expanduser('~'))
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
    text_content += [
        (N, N-1, 'pmc_oa', 'text', 'fulltext',
         zip_string('MEK phosphorylates ERK. EGFR activates SHC.'))
        ]
    db.copy('text_content', text_content,
            cols=('id', 'text_ref_id', 'source', 'format', 'text_type',
                  'content'))

    # Put an id file on s3
    basename = 'local_db_test_run'
    s3_prefix = 'reading_results/%s/' % basename
    s3.put_object(Bucket='bigmech', Key=s3_prefix + 'id_list',
                  Body='\n'.join(['tcid: %d' % i
                                  for i in range(len(text_content))]))

    # Call the reading tool
    sub = srp.DbReadingSubmitter(basename, ['sparser'])
    job_name, cmd = sub._make_command(0, len(text_content))
    cmd += ['--test']
    check_call(cmd)
    sub.produce_report()

    # Remove garbage on s3
    res = s3.list_objects(Bucket='bigmech', Prefix=s3_prefix)
    for entry in res['Contents']:
        print("Removing %s..." % entry['Key'])
        s3.delete_object(Bucket='bigmech', Key=entry['Key'])
    return


@attr('nonpublic')
def test_normal_pmid_reading_call():
    chdir(path.expanduser('~'))
    # Put an id file on s3
    basename = 'local_pmid_test_run'
    s3_prefix = 'reading_results/%s/' % basename
    s3.put_object(Bucket='bigmech', Key=s3_prefix + 'pmids',
                  Body='\n'.join(['PMID000test%d' % n for n in range(4)]))

    # Call the reading tool
    sub = srp.PmidSubmitter(basename, ['sparser'])
    job_name, cmd = sub._make_command(0, 2)
    check_call(cmd)

    # Remove garbage on s3
    res = s3.list_objects(Bucket='bigmech', Prefix=s3_prefix)
    for entry in res['Contents']:
        print("Removing %s..." % entry['Key'])
        s3.delete_object(Bucket='bigmech', Key=entry['Key'])
    return
