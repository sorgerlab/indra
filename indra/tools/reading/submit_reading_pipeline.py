from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import boto3
import logging
import botocore.session
from time import sleep
from datetime import datetime
from indra.literature import elsevier_client as ec
from indra.literature.elsevier_client import _ensure_api_keys
from indra.util.aws import get_job_log, tag_instance, get_batch_command

bucket_name = 'bigmech'

logger = logging.getLogger('aws_reading')


class BatchReadingError(Exception):
    pass


def wait_for_complete(queue_name, job_list=None, job_name_prefix=None,
                      poll_interval=10, idle_log_timeout=None,
                      kill_on_log_timeout=False, stash_log_method=None,
                      tag_instances=False):
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
    idle_log_timeout : Optional[int] or None
        If not None, then track the logs of the active jobs, and if new output
        is not produced after `idle_log_timeout` seconds, a warning is printed.
        If `kill_on_log_timeout` is set to True, the job will also be
        terminated.
    kill_on_log_timeout : Optional[bool]
        If True, and if `idle_log_timeout` is set, jobs will be terminated
        after timeout. This has no effect if `idle_log_timeout` is None.
        Default is False.
    stash_log_method : Optional[str]
        Select a method to store the job logs, either 's3' or 'local'. If no
        method is specified, the logs will not be loaded off of AWS. If 's3' is
        specified, then `job_name_prefix` must also be given, as this will
        indicate where on s3 to store the logs.
    """
    if stash_log_method == 's3' and job_name_prefix is None:
        raise Exception('A job_name_prefix is required to post logs on s3.')

    start_time = datetime.now()
    if job_list is None:
        job_id_list = []
    else:
        job_id_list = [job['jobId'] for job in job_list]

    def get_jobs_by_status(status, job_id_filter=None, job_name_prefix=None):
        res = batch_client.list_jobs(jobQueue=queue_name,
                                     jobStatus=status, maxResults=10000)
        jobs = res['jobSummaryList']
        if job_name_prefix:
            jobs = [job for job in jobs if
                    job['jobName'].startswith(job_name_prefix)]
        if job_id_filter:
            jobs = [job_def for job_def in jobs
                    if job_def['jobId'] in job_id_filter]
        return jobs

    job_log_dict = {}

    def check_logs(job_defs):
        stalled_jobs = set()
        for job_def in job_defs:
            try:
                log_lines = get_job_log(job_def, write_file=False)
                jid = job_def['jobId']
                now = datetime.now()
                if jid not in job_log_dict.keys():
                    logger.info("Adding job %s to the log tracker at %s."
                                % (jid, now))
                    job_log_dict[jid] = {'log': log_lines,
                                         'check_time': now}
                elif len(job_log_dict[jid]['log']) == len(log_lines):
                    check_dt = now - job_log_dict[jid]['check_time']
                    logger.warning(('Job \'%s\' has not produced output for '
                                    '%d seconds.')
                                   % (job_def['jobName'], check_dt.seconds))
                    if check_dt.seconds > idle_log_timeout:
                        logger.warning("Job \'%s\' has stalled."
                                       % job_def['jobName'])
                        stalled_jobs.add(jid)
                else:
                    old_log = job_log_dict[jid]['log']
                    old_log += log_lines[len(old_log):]
                    job_log_dict[jid]['check_time'] = now
            except Exception as e:
                logger.error("Failed to check log for: %s" % str(job_def))
                logger.exception(e)
        return stalled_jobs

    # Don't start watching jobs added after this command was initialized.
    observed_job_def_set = set()

    if stash_log_method is not None:
        def get_set_of_job_tuples(job_defs):
            return {tuple([(k, v) for k, v in job_def.items()
                           if k in ['jobName', 'jobId']])
                    for job_def in job_defs}

    batch_client = boto3.client('batch')
    if tag_instances:
        ecs_cluster_name = get_ecs_cluster_for_queue(queue_name, batch_client)

    terminate_msg = 'Job log has stalled for at least %f minutes.'
    terminated_jobs = set()
    while True:
        pre_run = []
        for status in ('SUBMITTED', 'PENDING', 'RUNNABLE', 'STARTING'):
            pre_run += get_jobs_by_status(status, job_id_list, job_name_prefix)
        running = get_jobs_by_status('RUNNING', job_id_list, job_name_prefix)
        failed = get_jobs_by_status('FAILED', job_id_list, job_name_prefix)
        done = get_jobs_by_status('SUCCEEDED', job_id_list, job_name_prefix)

        if stash_log_method is not None:
            observed_job_def_set |= get_set_of_job_tuples(pre_run + running)

        logger.info('(%d s)=(pre: %d, running: %d, failed: %d, done: %d)' %
                    ((datetime.now() - start_time).seconds, len(pre_run),
                     len(running), len(failed), len(done)))

        # Check the logs for new output, and possibly terminate some jobs.
        if idle_log_timeout is not None:
            stalled_jobs = check_logs(running)
            if kill_on_log_timeout:
                # Keep track of terminated jobs so we don't send a terminate
                # message twice.
                for jid in stalled_jobs.difference(terminated_jobs):
                    batch_client.terminate_job(
                        jobId=jid,
                        reason=terminate_msg % (idle_log_timeout/60.0)
                        )
                    logger.info('Terminating %s.' % jid)
                    terminated_jobs.add(jid)

        if job_id_list:
            if (len(failed) + len(done)) == len(job_id_list):
                ret = 0
                break
        else:
            if (len(failed) + len(done) > 0) and \
               (len(pre_run) + len(running) == 0):
                ret = 0
                break

        if tag_instances:
            tag_instances_on_cluster(ecs_cluster_name)
        sleep(poll_interval)

    # Stash the logs
    if stash_log_method:
        success_ids = [job_def['jobId'] for job_def in done]
        failure_ids = [job_def['jobId'] for job_def in failed]
        stash_logs(observed_job_def_set, success_ids, failure_ids, queue_name,
                   stash_log_method, job_name_prefix,
                   start_time.strftime('%Y%m%d_%H%M%S'))

    return ret


def stash_logs(job_defs, success_ids, failure_ids, queue_name, method='local',
               job_name_prefix=None, tag='stash'):
    if method == 's3':
        s3_client = boto3.client('s3')

        def stash_log(log_str, name_base):
            name = '%s_%s.log' % (name_base, tag)
            s3_client.put_object(
                Bucket=bucket_name,
                Key='reading_results/%s/logs/%s/%s' % (
                    job_name_prefix,
                    queue_name,
                    name),
                Body=log_str
                )

    elif method == 'local':
        if job_name_prefix is None:
            job_name_prefix = 'batch_%s' % tag
        dirname = '%s_job_logs' % job_name_prefix
        os.mkdir(dirname)

        def stash_log(log_str, name_base):
            with open(os.path.join(dirname, name_base + '.log'), 'w') as f:
                f.write(log_str)

    for job_def_tpl in job_defs:
        try:
            job_def = dict(job_def_tpl)
            lines = get_job_log(job_def, write_file=False)
            if lines is None:
                logger.warning("No logs found for %s." % job_def['jobName'])
                continue
            log_str = ''.join(lines)
            base_name = job_def['jobName']
            if job_def['jobId'] in success_ids:
                base_name += '_SUCCESS'
            elif job_def['jobId'] in failure_ids:
                base_name += '_FAILED'
            logger.info('Stashing ' + base_name)
            stash_log(log_str, base_name)
        except Exception as e:
            logger.error("Failed to save logs for: %s" % str(job_def_tpl))
            logger.exception(e)
    return


def get_ecs_cluster_for_queue(queue_name, batch_client=None):
    """Get the name of the ecs cluster using the batch client."""
    if batch_client is None:
        batch_client = boto3.client('batch')

    queue_resp = batch_client.describe_job_queues(jobQueues=[queue_name])
    if len(queue_resp['jobQueues']) == 1:
        queue = queue_resp['jobQueues'][0]
    else:
        raise BatchReadingError('Error finding queue with name %s.'
                                % queue_name)

    compute_env_names = queue['computeEnvironmentOrder']
    if len(compute_env_names) == 1:
        compute_env_name = compute_env_names[0]['computeEnvironment']
    else:
        raise BatchReadingError('Error finding the compute environment name '
                                'for %s.' % queue_name)

    compute_envs = batch_client.describe_compute_environments(
        computeEnvironments=[compute_env_name]
        )['computeEnvironments']
    if len(compute_envs) == 1:
        compute_env = compute_envs[0]
    else:
        raise BatchReadingError("Error getting compute environment %s for %s. "
                                "Got %d enviornments instead of 1."
                                % (compute_env_name, queue_name,
                                   len(compute_envs)))

    ecs_cluster_name = os.path.basename(compute_env['ecsClusterArn'])
    return ecs_cluster_name


def tag_instances_on_cluster(cluster_name, project='cwc'):
    """Adds project tag to untagged instances in a given cluster.

    Parameters
    ----------
    cluster_name : str
        The name of the AWS ECS cluster in which running instances
        should be tagged.
    project : str
        The name of the project to tag instances with.
    """
    # Get the relevent instance ids from the ecs cluster
    ecs = boto3.client('ecs')
    task_arns = ecs.list_tasks(cluster=cluster_name)['taskArns']
    if not task_arns:
        return
    tasks = ecs.describe_tasks(cluster=cluster_name, tasks=task_arns)['tasks']
    container_instances = ecs.describe_container_instances(
        cluster=cluster_name,
        containerInstances=[task['containerInstanceArn'] for task in tasks]
        )['containerInstances']
    ec2_instance_ids = [ci['ec2InstanceId'] for ci in container_instances]

    # Instantiate each instance to tag as a resource and create project tag
    for instance_id in ec2_instance_ids:
        tag_instance(instance_id, project=project)
    return


@_ensure_api_keys('remote batch reading', [])
def get_elsevier_api_keys():
    return [
        {'name': ec.API_KEY_ENV_NAME,
         'value': ec.ELSEVIER_KEYS.get('X-ELS-APIKey', '')},
        {'name': ec.INST_KEY_ENV_NAME,
         'value': ec.ELSEVIER_KEYS.get('X-ELS-Insttoken', '')},
        ]


def get_environment():
    # Get AWS credentials
    # http://stackoverflow.com/questions/36287720/boto3-get-credentials-dynamically
    session = botocore.session.get_session()
    access_key = session.get_credentials().access_key
    secret_key = session.get_credentials().secret_key

    # Get the Elsevier keys from the Elsevier client
    environment_vars = [
        {'name': 'AWS_ACCESS_KEY_ID',
         'value': access_key},
        {'name': 'AWS_SECRET_ACCESS_KEY',
         'value': secret_key}
        ]
    environment_vars += get_elsevier_api_keys()

    # Only include values that are not empty.
    return [var_dict for var_dict in environment_vars
            if var_dict['value'] and var_dict['name']]


def submit_reading(basename, pmid_list_filename, readers, start_ix=None,
                   end_ix=None, pmids_per_job=3000, num_tries=2,
                   force_read=False, force_fulltext=False, project_name=None):
    # Upload the pmid_list to Amazon S3
    pmid_list_key = 'reading_results/%s/pmids' % basename
    s3_client = boto3.client('s3')
    s3_client.upload_file(pmid_list_filename, bucket_name, pmid_list_key)

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
        command_list = get_batch_command(
            ['python', '-m', 'indra.tools.reading.pmid_reading.read_pmids_aws',
             basename, '/tmp', '16', str(job_start_ix), str(job_end_ix), '-r']
            + readers,
            purpose='pmid_reading',
            project=project_name
            )
        if force_read:
            command_list.append('--force_read')
        if force_fulltext:
            command_list.append('--force_fulltext')
        logger.info('Commands list: %s' % str(command_list))
        job_info = batch_client.submit_job(
            jobName=job_name,
            jobQueue='run_reach_queue',
            jobDefinition='run_reach_jobdef',
            containerOverrides={
                'environment': environment_vars,
                'command': command_list},
            retryStrategy={'attempts': num_tries}
            )
        logger.info("submitted...")
        job_list.append({'jobId': job_info['jobId']})
    return job_list


def submit_db_reading(basename, id_list_filename, readers, start_ix=None,
                      end_ix=None, pmids_per_job=3000, num_tries=2,
                      force_read=False, force_fulltext=False,
                      read_all_fulltext=False, project_name=None,
                      max_reach_input_len=None, max_reach_space_ratio=None):
    # Upload the pmid_list to Amazon S3
    pmid_list_key = 'reading_inputs/%s/id_list' % basename
    s3_client = boto3.client('s3')
    s3_client.upload_file(id_list_filename, bucket_name, pmid_list_key)

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

    if force_read:
        mode = 'all'
    else:
        mode = 'unread'

    # Iterate over the list of PMIDs and submit the job in chunks
    batch_client = boto3.client('batch', region_name='us-east-1')
    job_list = []
    for job_start_ix in range(start_ix, end_ix, pmids_per_job):
        job_end_ix = job_start_ix + pmids_per_job
        if job_end_ix > end_ix:
            job_end_ix = end_ix
        job_name = '%s_%d_%d' % (basename, job_start_ix, job_end_ix)
        command_list = get_batch_command(
            ['python', '-m', 'indra.tools.reading.db_reading.read_db_aws',
             basename, '/tmp', mode, '32', str(job_start_ix), str(job_end_ix),
             '-r'] + readers,
            purpose='db_reading',
            project=project_name
            )
        if force_fulltext:
            command_list.append('--force_fulltext')
        if read_all_fulltext:
            command_list.append('--read_all_fulltext')
        if max_reach_input_len is not None:
            command_list += ['--max_reach_input_len', max_reach_input_len]
        if max_reach_space_ratio is not None:
            command_list += ['--max_reach_space_ratio', max_reach_space_ratio]
        logger.info('Command list: %s' % str(command_list))
        job_info = batch_client.submit_job(
            jobName=job_name,
            jobQueue='run_db_reading_queue',
            jobDefinition='run_db_reading_jobdef',
            containerOverrides={
                'environment': environment_vars,
                'command': command_list},
            retryStrategy={'attempts': num_tries}
            )
        logger.info("submitted...")
        job_list.append({'jobId': job_info['jobId']})
    return job_list


def submit_combine(basename, readers, job_ids=None, project_name=None):
    if job_ids is not None and len(job_ids) > 20:
        print("WARNING: boto3 cannot support waiting for more than 20 jobs.")
        print("Please wait for the reading to finish, then run again with the")
        print("`combine` option.")
        return

    # Get environment variables
    environment_vars = get_environment()

    job_name = '%s_combine_reading_results' % basename
    command_list = get_batch_command(
        ['python', '-m', 'indra.tools.reading.assemble_reading_stmts_aws',
         basename, '-r'] + readers,
        purpose='pmid_reading',
        project=project_name
        )
    logger.info('Command list: %s' % str(command_list))
    kwargs = {'jobName': job_name, 'jobQueue': 'run_reach_queue',
              'jobDefinition': 'run_reach_jobdef',
              'containerOverrides': {'environment': environment_vars,
                                     'command': command_list,
                                     'memory': 60000, 'vcpus': 1}}
    if job_ids:
        kwargs['dependsOn'] = job_ids
    batch_client = boto3.client('batch')
    batch_client.submit_job(**kwargs)
    logger.info("submitted...")


if __name__ == '__main__':
    import argparse

    # Create the top-level parser
    parser = argparse.ArgumentParser(
        'submit_reading_pipeline_aws.py',
        description=('Run reading with either the db or remote resources. For '
                     'more specific help, select one of the Methods with the '
                     '`-h` option.'),
        epilog=('Note that `python wait_for_complete.py ...` should be run as '
                'soon as this command completes successfully. For more '
                'details use `python wait_for_complete.py -h`.')
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
        nargs='+',
        help='Choose which reader(s) to use.'
        )
    parent_submit_parser.add_argument(
        '--project',
        help=('Set the project name. Default is DEFAULT_AWS_PROJECT in the '
              'config.')
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
    parent_db_parser.add_argument(
        '--read_all_fulltext',
        action='store_true',
        help='Read all available fulltext for input ids, not just the best.'
        )
    parent_db_parser.add_argumet(
        '--max_reach_space_ratio',
        type=float,
        help='Set the maximum ratio of spaces to non-spaces for REACH input.',
        default=None
    )
    parent_db_parser.add_argument(
        '--max_reach_input_len',
        type=int,
        help='Set the maximum length of content that REACH will read.',
        default=None
    )

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

    job_ids = None
    if args.method == 'no-db':
        if args.job_type in ['read', 'full']:
            job_ids = submit_reading(
                args.basename,
                args.input_file,
                args.readers,
                args.start_ix,
                args.end_ix,
                args.ids_per_job,
                2,
                args.force_read,
                args.force_fulltext,
                args.project
                )
        if args.job_type in ['combine', 'full']:
            submit_combine(args.basename, args.readers, job_ids, args.project)
    elif args.method == 'with-db':
        job_ids = submit_db_reading(
            args.basename,
            args.input_file,
            args.readers,
            args.start_ix,
            args.end_ix,
            args.ids_per_job,
            2,
            args.force_read,
            args.force_fulltext,
            args.read_all_fulltext,
            args.project,
            args.max_reach_input_len,
            args.max_reach_space_ratio
            )
