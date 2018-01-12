from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import boto3
import logging
import botocore.session
from time import sleep
from indra.literature import elsevier_client as ec
from indra.tools.reading.read_pmids import READER_DICT
from docutils.nodes import description

bucket_name = 'bigmech'

logger = logging.getLogger('aws_reading')


def wait_for_complete(queue_name, job_list=None, job_name_prefix=None,
                      poll_interval=10):
    """Return when all jobs in the given list finished.

    If not job list is given, return when all jobs in queue finished.

    Parameters
    ----------
    queue_name : str
        The name of the queue to wait for completion.
    job_list : Optional[list(dict)]
        A list of jobID-s in a dict, as returned by the submit function.
        Example: [{'jobId': 'e6b00f24-a466-4a72-b735-d205e29117b4'}, ...]
        If not given, this function will return if all jobs completed.
    job_name_prefix : Optional[str]
        A prefix for the name of the jobs to wait for. This is useful if the
        explicit job list is not available but filtering is needed.
    poll_interval : Optional[int]
        The time delay between API calls to check the job statuses.
    """
    if job_list is None:
        job_list = []

    def get_jobs_by_status(status, job_filter=None, job_name_prefix=None):
        res = batch_client.list_jobs(jobQueue=queue_name,
                                     jobStatus=status, maxResults=10000)
        jobs = res['jobSummaryList']
        if job_name_prefix:
            jobs = [job for job in jobs if
                    job['jobName'].startswith(job_name_prefix)]
        ids = [job['jobId'] for job in jobs]
        if job_filter:
            ids = [job_id for job_id in ids if job_id in job_filter]
        return ids

    job_list = [job['jobId'] for job in job_list]

    batch_client = boto3.client('batch')

    total_time = 0
    while True:
        pre_run = []
        for status in ('SUBMITTED', 'PENDING', 'RUNNABLE', 'STARTING'):
            pre_run += get_jobs_by_status(status, job_list, job_name_prefix)
        running = get_jobs_by_status('RUNNING', job_list, job_name_prefix)
        failed = get_jobs_by_status('FAILED', job_list, job_name_prefix)
        done = get_jobs_by_status('SUCCEEDED', job_list, job_name_prefix)

        logger.info('(%d s)=(pre: %d, running: %d, failed: %d, done: %d)' %
                    (total_time, len(pre_run), len(running),
                     len(failed), len(done)))

        if job_list:
            if (len(failed) + len(done)) == len(job_list):
                return 0
        else:
            if (len(failed) + len(done) > 0) and \
                (len(pre_run) + len(running) == 0):
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


def submit_db_reading(basename, id_list_filename, readers, start_ix=None,
                      end_ix=None, pmids_per_job=3000, num_tries=2):
    # Upload the pmid_list to Amazon S3
    pmid_list_key = 'reading_inputs/%s/id_list' % basename
    s3_client = boto3.client('s3')
    s3_client.upload_file(id_list_filename, 'bigmech', pmid_list_key)

    # If no end index is specified, read all the PMIDs
    if end_ix is None:
        with open(id_list_filename, 'rt') as f:
            lines = f.readlines()
            end_ix = len(lines)

    if start_ix is None:
        start_ix = 0

    # Get environment variables
    environment_vars = get_environment()

    # Fix reader options
    if 'all' in readers:
        readers = ['reach', 'sparser']

    # Iterate over the list of PMIDs and submit the job in chunks
    batch_client = boto3.client('batch')
    job_list = []
    for job_start_ix in range(start_ix, end_ix, pmids_per_job):
        job_end_ix = job_start_ix + pmids_per_job
        if job_end_ix > end_ix:
            job_end_ix = end_ix
        job_name = '%s_%d_%d' % (basename, job_start_ix, job_end_ix)
        command_list = ['python', '-m',
                        'indra.tools.reading.read_db_aws',
                        basename, '/tmp', 'unread', '32', str(job_start_ix),
                        str(job_end_ix), '-r'] + readers
        print(command_list)
        job_info = batch_client.submit_job(
            jobName=job_name,
            jobQueue='run_db_reading_queue',
            jobDefinition='run_db_reading_jobdef',
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
        description='Run reading with either the db or remote resources.'
        )
    subparsers = parser.add_subparsers(title='Method')
    subparsers.required = True
    subparsers.dest = 'method'

    # Create parser class for first layer of options
    grandparent_reading_parser = argparse.ArgumentParser(
        description='Run machine reading using AWS Batch.',
        add_help=False
        )

    # Create parent parser classes for second layer of options
    parent_submit_parser = argparse.ArgumentParser(add_help=False)
    parent_submit_parser.add_argument(
        'basename',
        help='Defines job names and S3 keys'
        )
    parent_submit_parser.add_argument(
        '-r', '--readers',
        dest='readers',
        choices=['sparser', 'reach', 'all'],
        default=['all'],
        help='Choose which reader(s) to use.'
        )
    parent_read_parser = argparse.ArgumentParser(add_help=False)
    parent_read_parser.add_argument(
        'input_file',
        help=('Path to file containing input ids of content to read. For the '
              'no-db options, this is simply a file with each line being a '
              'pmid. For the with-db options, this is a file where each line '
              'is of the form \'<id type>:<id>\', for example \'pmid:12345\'')
        )
    parent_read_parser.add_argument(
        '--start_ix',
        type=int,
        help='Start index of ids to read.'
        )
    parent_read_parser.add_argument(
        '--end_ix',
        type=int,
        help='End index of ids to read. If `None`, read content from all ids.'
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
        '--ids_per_job',
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
    parent_db_parser = argparse.ArgumentParser(add_help=False)
    '''Not currently supported
    parent_db_parser.add_argument(
        '--no_upload',
        action='store_true',
        help='Don\'t upload results to the database.'
        )
    '''

    # Make non_db_parser and get subparsers
    non_db_parser = subparsers.add_parser(
        'no-db',
        parents=[grandparent_reading_parser],
        description=('Run reading by collecting content, and save as pickles. '
                     'This option requires that ids are given as a list of '
                     'pmids, one line per pmid.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    non_db_subparsers = non_db_parser.add_subparsers(
        title='Job Type',
        help='Type of jobs to submit.'
        )
    non_db_subparsers.required = True
    non_db_subparsers.dest = 'job_type'

    # Create subparsers for the no-db option.
    read_parser = non_db_subparsers.add_parser(
        'read',
        parents=[parent_read_parser, parent_submit_parser],
        help='Run REACH and cache INDRA Statements on S3.',
        description='Run REACH and cache INDRA Statements on S3.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    combine_parser = non_db_subparsers.add_parser(
        'combine',
        parents=[parent_submit_parser],
        help='Combine INDRA Statement subsets into a single file.',
        description='Combine INDRA Statement subsets into a single file.'
        )
    full_parser = non_db_subparsers.add_parser(
        'full',
        parents=[parent_read_parser, parent_submit_parser],
        help='Run REACH and combine INDRA Statements when done.',
        description='Run REACH and combine INDRA Statements when done.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    # Make db parser and get subparsers.
    db_parser = subparsers.add_parser(
        'with-db',
        parents=[grandparent_reading_parser, parent_submit_parser,
                 parent_read_parser, parent_db_parser],
        description=('Run reading with content on the db and submit results. '
                     'In this option, ids in \'input_file\' are given in the '
                     'format \'<id type>:<id>\'. Unlike no-db, there is no '
                     'need to combine pickles, and therefore no need to '
                     'specify your task further.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    args = parser.parse_args()

    if not isinstance(args.readers, (list, tuple)):
        args.readers = [args.readers]

    job_ids = None
    if args.method == 'no-db':
        if args.job_type in ['read', 'full']:
            job_ids = submit_reading(
                args.basename,
                args.input_file,
                args.readers,
                args.start_ix,
                args.end_ix,
                args.ids_per_job
                )
        if args.job_type in ['combine', 'full']:
            submit_combine(args.basename, args.readers, job_ids)
    elif args.method == 'with-db':
        job_ids = submit_db_reading(
            args.basename,
            args.input_file,
            args.readers,
            args.start_ix,
            args.end_ix,
            args.ids_per_job
            )
