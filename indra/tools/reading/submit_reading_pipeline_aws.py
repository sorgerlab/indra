from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import boto3
import logging
import botocore.session
from time import sleep
from indra.literature import elsevier_client as ec
from indra.tools.reading.read_pmids import READER_DICT

bucket_name = 'bigmech'

logger = logging.getLogger('aws_reading')


def wait_for_complete(job_list, poll_interval=10):
    """Return when all jobs in the given list finished.

    If not job list is given, return when all jobs in queue finished.

    Parameters
    ----------
    job_list : Optional[list(dict)]
        A list of jobID-s in a dict, as returned by the submit function.
        Example: [{'jobId': 'e6b00f24-a466-4a72-b735-d205e29117b4'}, ...]
        If not given, this function will return if all jobs completed.
    poll_interval : Optional[int]
        The time delay between API calls to check the job statuses.
    """
    def get_jobs_by_status(status, job_filter=None):
        res = batch_client.list_jobs(jobQueue='run_reach_queue',
                                     jobStatus=status)
        ids = [job['jobId'] for job in res['jobSummaryList']]
        if job_filter:
            ids = [job_id for job_id in ids if job_id in job_filter]
        return ids

    job_list = [job['jobId'] for job in job_list]

    batch_client = boto3.client('batch')

    total_time = 0
    while True:
        not_done = []
        for status in ('SUBMITTED', 'PENDING', 'RUNNABLE', 'STARTING',
                       'RUNNING'):
            not_done += get_jobs_by_status(status, job_list)
        failed = get_jobs_by_status('FAILED', job_list)
        done = get_jobs_by_status('SUCCEEDED', job_list)

        logger.info(
            '(%d s)=(not done: %d, failed: %d, done: %d)' %
            (total_time, len(not_done), len(failed), len(done))
            )

        if job_list:
            if (len(failed) + len(done)) == len(job_list):
                return 0
        else:
            if (len(failed) + len(done) > 0) and (len(not_done) == 0):
                return 0

        tag_instances()
        sleep(poll_interval)
        total_time += poll_interval


def tag_instances(project='bigmechanism'):
    """Adds project tag to untagged fleet instances."""
    # First, get all the instances
    ec2_client = boto3.client('ec2')
    resp = ec2_client.describe_instances()
    instances = []
    for res in resp.get('Reservations', []):
        instances += res.get('Instances', [])
    instances_to_tag = []
    # Check each instance to see if it's tagged and if it's a spot fleet
    # instance
    for instance in instances:
        tagged = False
        need_tag = False
        for tag in instance.get('Tags', []):
            if tag.get('Key') == 'project':
                tagged = True
            elif tag.get('Key') == 'aws:ec2spot:fleet-request-id':
                need_tag = True
        if not tagged and need_tag:
            instances_to_tag.append(instance['InstanceId'])
    # Instantiate each instance to tag as a resource and create project tag
    ec2 = boto3.resource('ec2')
    for instance_id in instances_to_tag:
        logger.info('Adding project tag to instance %s' % instance_id)
        instance = ec2.Instance(instance_id)
        instance.create_tags(Tags=[{'Key': 'project',
                                    'Value': project}])


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
         'value': secret_key}
        ]
    return environment_vars


def submit_reading(basename, pmid_list_filename, readers, start_ix=None,
                   end_ix=None, pmids_per_job=3000, num_tries=2):
    # Upload the pmid_list to Amazon S3
    pmid_list_key = 'reading_results/%s/pmids' % basename
    s3_client = boto3.client('s3')
    s3_client.upload_file(pmid_list_filename, 'bigmech', pmid_list_key)

    # If no end index is specified, read all the PMIDs
    if end_ix is None:
        with open(pmid_list_filename, 'rt') as f:
            lines = f.readlines()
            end_ix = len(lines)

    if start_ix is None:
        start_ix = 0

    # Get environment variables
    environment_vars = get_environment()

    # Iterate over the list of PMIDs and submit the job in chunks
    batch_client = boto3.client('batch')
    job_list = []
    for job_start_ix in range(start_ix, end_ix, pmids_per_job):
        job_end_ix = job_start_ix + pmids_per_job
        if job_end_ix > end_ix:
            job_end_ix = end_ix
        job_name = '%s_%d_%d' % (basename, job_start_ix, job_end_ix)
        command_list = ['python', '-m',
                        'indra.tools.reading.read_pmids_aws',
                        basename, '/tmp', '16', str(job_start_ix),
                        str(job_end_ix), '-r'] + readers
        print(command_list)
        job_info = batch_client.submit_job(
            jobName=job_name,
            jobQueue='run_reach_queue',
            jobDefinition='run_reach_jobdef',
            containerOverrides={
                'environment': environment_vars,
                'command': command_list},
            retryStrategy={'attempts': num_tries}
            )
        print("submitted...")
        job_list.append({'jobId': job_info['jobId']})
    return job_list


def submit_combine(basename, readers, job_ids=None):
    if job_ids is not None and len(job_ids) > 20:
        print("WARNING: boto3 cannot support waiting for more than 20 jobs.")
        print("Please wait for the reading to finish, then run again with the")
        print("`combine` option.")
        return

    # Get environment variables
    environment_vars = get_environment()

    job_name = '%s_combine_reading_results' % basename
    command_list = ['python', '-m',
                    'indra.tools.reading.assemble_reading_stmts_aws',
                    basename, '-r'] + readers
    print(command_list)
    kwargs = {'jobName': job_name, 'jobQueue': 'run_reach_queue',
              'jobDefinition': 'run_reach_jobdef',
              'containerOverrides': {'environment': environment_vars,
                                     'command': command_list,
                                     'memory': 60000, 'vcpus': 1}}
    if job_ids:
        kwargs['dependsOn'] = job_ids
    batch_client = boto3.client('batch')
    batch_client.submit_job(**kwargs)
    print("submitted...")


if __name__ == '__main__':
    import argparse

    # Create the top-level parser
    parser = argparse.ArgumentParser(
        'submit_reading_pipeline_aws.py',
        description='Run machine reading with REACH using AWS Batch.'
        )
    subparsers = parser.add_subparsers(
        title='Job Type',
        help='Type of jobs to submit.'
        )
    subparsers.required = True
    subparsers.dest = 'job_type'
    parent_submit_parser = argparse.ArgumentParser(add_help=False)
    parent_submit_parser.add_argument(
        'basename',
        help='Defines job names and S3 keys'
        )
    parent_submit_parser.add_argument(
        '-r', '--readers',
        dest='readers',
        choices=list(READER_DICT.keys()) + ['all'],
        default=['all'],
        nargs='+',
        help='Choose which reader(s) to use.'
        )
    parent_read_parser = argparse.ArgumentParser(add_help=False)
    parent_read_parser.add_argument(
        'pmid_file',
        help='Path to file containing PMIDs to read'
        )
    parent_read_parser.add_argument(
        '--start_ix',
        type=int,
        help='Start index of PMIDs to read.'
        )
    parent_read_parser.add_argument(
        '--end_ix',
        type=int,
        help='End index of PMIDs to read (default: read all PMIDs)'
        )
    parent_read_parser.add_argument(
        '--force_read',
        action='store_true',
        help='Read papers even if previously read by current REACH.'
        )
    parent_read_parser.add_argument(
        '--force_fulltext',
        action='store_true',
        help='Get full text content even if content already on S3.'
        )
    parent_read_parser.add_argument(
        '--pmids_per_job',
        default=3000,
        type=int,
        help='Number of PMIDs to read for each AWS Batch job.'
        )
    ''' Not currently supported.
    parent_read_parser.add_argument(
        '--num_tries',
        default=2,
        type=int,
        help='Maximum number of times to try running job.'
        )
    '''
    read_parser = subparsers.add_parser(
        'read',
        parents=[parent_read_parser, parent_submit_parser],
        help='Run REACH and cache INDRA Statements on S3.',
        description='Run REACH and cache INDRA Statements on S3.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    combine_parser = subparsers.add_parser(
        'combine',
        parents=[parent_submit_parser],
        help='Combine INDRA Statement subsets into a single file.',
        description='Combine INDRA Statement subsets into a single file.'
        )
    full_parser = subparsers.add_parser(
        'full',
        parents=[parent_read_parser, parent_submit_parser],
        help='Run REACH and combine INDRA Statements when done.',
        description='Run REACH and combine INDRA Statements when done.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    args = parser.parse_args()

    job_ids = None
    if args.job_type in ['read', 'full']:
        job_ids = submit_reading(
            args.basename,
            args.pmid_file,
            args.readers,
            args.start_ix,
            args.end_ix,
            args.pmids_per_job
            )
    if args.job_type in ['combine', 'full']:
        submit_combine(args.basename, args.readers, job_ids)
