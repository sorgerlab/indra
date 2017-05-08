from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

def get_environment():
    # Get AWS credentials
    # http://stackoverflow.com/questions/36287720/boto3-get-credentials-dynamically
    session = botocore.session.get_session()
    access_key = session.get_credentials().access_key
    secret_key = session.get_credentials().secret_key

    # Get the Elsevier keys from the Elsevier client
    environment_vars = [
            {'name': ec.api_key_env_name,
             'value': ec.elsevier_keys.get('X-ELS-APIKey')},
            {'name': ec.inst_key_env_name,
             'value': ec.elsevier_keys.get('X-ELS-Insttoken')},
            {'name': 'AWS_ACCESS_KEY_ID',
             'value': access_key},
            {'name': 'AWS_SECRET_ACCESS_KEY',
             'value': secret_key}]
    return environment_vars


def submit_run_reach(basename, pmid_list_filename, start_ix=0, end_ix=None,
                 pmids_per_job=3000):
    # Upload the pmid_list to Amazon S3
    pmid_list_key = 'reading_results/%s/pmids' % basename
    s3_client = boto3.client('s3')
    s3_client.upload_file(pmid_list_filename, 'bigmech', pmid_list_key)

    # If no end index is specified, read all the PMIDs
    if end_ix is None:
        with open(pmid_list_filename, 'rt') as f:
            lines = f.readlines()
            end_ix = len(lines)

    # Get environment variables
    environment_vars = get_environment()

    # Iterate over the list of PMIDs and submit the job in chunks
    job_list = []
    for job_start_ix in range(start_ix, end_ix, pmids_per_job):
        job_end_ix = job_start_ix + pmids_per_job
        if job_end_ix > end_ix:
            job_end_ix = end_ix
        job_name = '%s_%d_%d' % (basename, job_start_ix, job_end_ix)
        command_list = ['python', '-m',
                        'indra.tools.reading.run_reach_on_pmids_aws',
                        basename, '/tmp', '16', str(job_start_ix),
                        str(job_end_ix)]
        print(command_list)
        job_info = batch_client.submit_job(jobName=job_name,
                       jobQueue='run_reach_queue',
                       jobDefinition='run_reach_jobdef',
                       containerOverrides={'environment': environment_vars,
                                           'command': command_list})
        job_list.append({'jobId': job_info['jobId']})
    return job_list

def submit_combine(basename, job_ids=None):
    # Get environment variables
    environment_vars = get_environment()

    job_name = '%s_combine_reach' % basename
    command_list = ['python', '-m',
                    'indra.tools.reading.assemble_reach_stmts_aws',
                    basename]
    print(command_list)
    kwargs = {'jobName': job_name, 'jobQueue': 'run_reach_queue',
              'jobDefinition': 'run_reach_jobdef',
              'containerOverrides': {'environment': environment_vars,
                                     'command': command_list,
                                     'memory': 60000, 'vcpus': 1}}
    if job_ids:
        kwargs['dependsOn'] = job_ids
    batch_client.submit_job(**kwargs)

if __name__ == '__main__':
    import sys
    import boto3
    import botocore.session
    from indra.literature import elsevier_client as ec

    bucket_name = 'bigmech'
    job_type =sys.argv[1]

    # Submit the reading job
    batch_client = boto3.client('batch')

    if job_type == 'read':
        basename = sys.argv[2]
        pmid_list_filename = sys.argv[3]
        start_ix = int(sys.argv[4])
        end_ix = int(sys.argv[5])
        pmids_per_job = int(sys.argv[6])
        # force_read
        # force_fulltext
        job_ids = submit_run_reach(basename, pmid_list_filename, start_ix,
                               end_ix, pmids_per_job)
    elif job_type == 'combine':
        basename = sys.argv[2]
        submit_combine(basename)
    elif job_type == 'full':
        basename = sys.argv[2]
        pmid_list_filename = sys.argv[3]
        start_ix = int(sys.argv[4])
        end_ix = int(sys.argv[5])
        pmids_per_job = int(sys.argv[6])
        job_ids = submit_run_reach(basename, pmid_list_filename, start_ix,
                                   end_ix, pmids_per_job)
        submit_combine(basename, job_ids)
    else:
        print('job type must be "read", "combine", or "full"')
        sys.exit(1)


