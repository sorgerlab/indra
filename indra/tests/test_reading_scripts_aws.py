
import boto3
from os import path, chdir
from subprocess import check_call
from nose.plugins.attrib import attr

from indra.tools.reading import submit_reading_pipeline as srp
from indra.sources import sparser


s3 = boto3.client('s3')


HERE = path.dirname(path.abspath(__file__))


@attr('nonpublic', 'notravis')
def test_normal_pmid_reading_call():
    chdir(path.expanduser('~'))
    # Put an id file on s3
    basename = 'local_pmid_test_run'
    s3_prefix = 'reading_results/%s/' % basename
    s3.put_object(Bucket='bigmech', Key=s3_prefix + 'pmids',
                  Body='\n'.join(['PMID000test%d' % n for n in range(4)]))

    # Call the reading tool
    sub = srp.PmidSubmitter(basename, ['sparser'])
    for job_name, cmd, job_def in sub._iter_commands(0, 2):
        check_call(cmd)

    # Remove garbage on s3
    res = s3.list_objects(Bucket='bigmech', Prefix=s3_prefix)
    for entry in res['Contents']:
        print("Removing %s..." % entry['Key'])
        s3.delete_object(Bucket='bigmech', Key=entry['Key'])
    return


@attr('nonpublic', 'notravis')
def test_bad_sparser():
    txt = ('Disruption of the AP-1 binding site reversed the transcriptional '
           'responses seen with Fos and Jun.')
    sp = sparser.process_text(txt, timeout=1)
    assert sp is None, "Reading succeeded unexpectedly."
