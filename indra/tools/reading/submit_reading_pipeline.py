import os
import boto3
import logging
import argparse

from time import sleep
from threading import Thread
from datetime import datetime
from indra.literature import elsevier_client as ec
from indra.literature.elsevier_client import _ensure_api_keys
from indra.tools.reading.readers import get_reader_classes
from indra.util.aws import tag_instance, get_batch_command, kill_all, get_ids,\
    JobLog

bucket_name = 'bigmech'

logger = logging.getLogger('indra.tools.reading.submit_reading_pipeline')


class BatchReadingError(Exception):
    pass


class BatchMonitor(object):
    """A monitor for batch jobs.

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
    result_record : dict
        A dict which will be modified in place to record the results of the job
    """
    def __init__(self, queue_name, job_list=None, job_name_prefix=None,
                 result_record=None):

        self.start_time = datetime.now()
        self.queue_name = queue_name
        self.job_name_prefix = job_name_prefix
        self.job_list = job_list

        self.result_record = {} if result_record is None else result_record
        self.job_log_dict = {}

        # Don't start watching jobs added after this command was initialized.
        self.observed_job_def_dict = {}

        self.batch_client = boto3.client('batch')

        self.job_id_list = None
        return

    def watch_and_wait(self, poll_interval=10, idle_log_timeout=None,
                       kill_on_log_timeout=False, stash_log_method=None,
                       tag_instances=False, wait_for_first_job=False):
        """Return when all jobs are finished.

        If no job list was given, return when all jobs in queue finished.

        Parameters
        ----------
        poll_interval : Optional[int]
            The time delay between API calls to check the job statuses.
        idle_log_timeout : Optional[int] or None
            If not None, then track the logs of the active jobs, and if new
            output is not produced after `idle_log_timeout` seconds, a warning
            is printed. If `kill_on_log_timeout` is set to True, the job will
            also be terminated.
        kill_on_log_timeout : Optional[bool]
            If True, and if `idle_log_timeout` is set, jobs will be terminated
            after timeout. This has no effect if `idle_log_timeout` is None.
            Default is False.
        stash_log_method : Optional[str]
            Select a method to store the job logs, either 's3' or 'local'. If
            no method is specified, the logs will not be loaded off of AWS. If
            's3' is specified, then `job_name_prefix` must also be given, as
            this will indicate where on s3 to store the logs.
        tag_instances : bool
            Default is False. If True, apply tags to the instances. This is
            today typically done by each job, so in most cases this should not
            be needed.
        wait_for_first_job : bool
            Don't exit until at least one job has been found. This is good if
            you are monitoring jobs that are submitted periodically, but can be
            a problem if there is a chance you might call this when no jobs
            will ever be run.
        """
        logger.info("Given %s jobs to track"
                    % ('no' if self.job_list is None else len(self.job_list)))
        if stash_log_method == 's3' and self.job_name_prefix is None:
            raise Exception('A job_name_prefix is required to post logs on s3.')
        if tag_instances:
            ecs_cluster_name = \
                get_ecs_cluster_for_queue(self.queue_name, self.batch_client)
        else:
            ecs_cluster_name = None
        terminate_msg = 'Job log has stalled for at least %f minutes.'
        terminated_jobs = set()
        stashed_id_set = set()
        found_a_job = False
        while True:
            pre_run = []
            self.job_id_list = get_ids(self.job_list)
            logger.info("Specifically tracked jobs: %s"
                        % len(self.job_id_list))
            for status in ('SUBMITTED', 'PENDING', 'RUNNABLE', 'STARTING'):
                pre_run += self.get_jobs_by_status(status)
            running = self.get_jobs_by_status('RUNNING')
            failed = self.get_jobs_by_status('FAILED')
            done = self.get_jobs_by_status('SUCCEEDED')

            if len(pre_run + running):
                found_a_job = True

            self.observed_job_def_dict.update(
                self.get_dict_of_job_tuples(pre_run + running)
            )

            logger.info('(%d s)=(pre: %d, running: %d, failed: %d, done: %d)' %
                        ((datetime.now() - self.start_time).seconds,
                         len(pre_run), len(running), len(failed), len(done)))

            # Check the logs for new output, and possibly terminate some jobs.
            stalled_jobs = self.check_logs(running)
            if idle_log_timeout is not None:
                if kill_on_log_timeout:
                    # Keep track of terminated jobs so we don't send a
                    # terminate message twice.
                    for jid in stalled_jobs - terminated_jobs:
                        self.batch_client.terminate_job(
                            jobId=jid,
                            reason=terminate_msg % (idle_log_timeout/60.0)
                        )
                        logger.info('Terminating %s.' % jid)
                        terminated_jobs.add(jid)

            # Check for end-conditions.
            if found_a_job or not wait_for_first_job:
                if self.job_id_list:
                    if (len(failed) + len(done)) == len(self.job_id_list):
                        logger.info("Total failed and done equals number of "
                                    "original tracked jobs. Ending.")
                        ret = 0
                        break
                else:
                    if (len(failed) + len(done) > 0) and \
                            (len(pre_run) + len(running) == 0):
                        logger.info("No job_id_list, but there are new "
                                    "finished jobs and no running or pre-"
                                    "running jobs. Ending.")
                        ret = 0
                        break

            if tag_instances:
                tag_instances_on_cluster(ecs_cluster_name)

            # Stash the logs of things that have finished so far. Note that
            # jobs terminated in this round will not be picked up until the
            # next round.
            self.stash_logs(stash_log_method, done, failed, stashed_id_set)

            sleep(poll_interval)

        # Pick up any stragglers
        self.stash_logs(stash_log_method, done, failed, stashed_id_set)

        self.result_record['terminated'] = terminated_jobs
        self.result_record['failed'] = failed
        self.result_record['succeeded'] = done

        return ret

    def stash_logs(self, method, done, failed, stashed_ids):
        if not method:
            return
        stash_logs(self.observed_job_def_dict, done, failed, self.queue_name,
                   method, self.job_name_prefix,
                   self.start_time.strftime('%Y%m%d_%H%M%S'),
                   ids_stashed=stashed_ids)

    def get_jobs_by_status(self, status):
        res = self.batch_client.list_jobs(jobQueue=self.queue_name,
                                          jobStatus=status, maxResults=10000)
        jobs = res['jobSummaryList']
        if self.job_name_prefix:
            jobs = [job for job in jobs if
                    job['jobName'].startswith(self.job_name_prefix)]
        if self.job_id_list is not None:
            jobs = [job_def for job_def in jobs
                    if job_def['jobId'] in self.job_id_list]
        return jobs

    def check_logs(self, job_defs):
        """Updates the job_log_dict."""
        stalled_jobs = set()

        # Check the status of all the jobs we're tracking.
        for job_def in job_defs:
            try:
                # Get the job id.
                jid = job_def['jobId']
                now = datetime.utcnow()
                if jid not in self.job_log_dict.keys():
                    # If the job is new...
                    logger.info("Adding job %s to the log job_log at %s."
                                % (jid, now))
                    # Instantiate a new job_log.
                    job_log = JobLog(job_def)
                    self.job_log_dict[jid] = job_log
                else:
                    job_log = self.job_log_dict[jid]

                pre_len = len(job_log)
                job_log.get_lines()
                post_len = len(job_log)

                if pre_len == post_len:
                    # If the job log hasn't changed, announce as such, and
                    # check to see if it has been the same for longer than
                    # stall time.
                    check_dt = now - job_log.latest_timestamp
                    logger.warning(('Job \'%s\' has not produced output for '
                                    '%d seconds.')
                                   % (job_def['jobName'], check_dt.seconds))
                    if check_dt.seconds > self.idle_log_timeout:
                        logger.warning("Job \'%s\' has stalled."
                                       % job_def['jobName'])
                        stalled_jobs.add(jid)

                # Dump the log lines as we go, reducing RAM usage.
                job_log.dump()
                job_log.clear_lines()
            except Exception as e:
                # Sometimes due to sync et al. issues, a part of this will fail
                # Such things are usually transitory issues so we keep trying.
                logger.error("Failed to check log for: %s" % str(job_def))
                logger.exception(e)

        # Pass up the set of job id's for stalled jobs.
        return stalled_jobs

    @staticmethod
    def get_dict_of_job_tuples(job_defs):
        return {jdef['jobId']: [(k, jdef[k]) for k in ['jobName', 'jobId']]
                for jdef in job_defs}


def _get_job_ids_to_stash(job_def_list, stashed_id_set):
    return [job_def['jobId'] for job_def in job_def_list
            if job_def['jobId'] not in stashed_id_set]


def stash_logs(job_defs, success_jobs, failure_jobs, queue_name, method='local',
               job_name_prefix=None, tag='stash', ids_stashed=None):
    if ids_stashed is None:
        ids_stashed = set()

    success_ids = _get_job_ids_to_stash(success_jobs, ids_stashed)
    failure_ids = _get_job_ids_to_stash(failure_jobs, ids_stashed)
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
    else:
        raise ValueError('Invalid method: %s' % method)

    for jobId, job_def_tpl in job_defs.items():
        if jobId not in success_ids and jobId not in failure_ids:
            continue  # Logs aren't done and ready to be loaded.
        try:
            job_def = dict(job_def_tpl)
            lines = get_job_log(job_def, write_file=False)
            if lines is None:
                logger.warning("No logs found for %s." % job_def['jobName'])
                continue
            log_str = ''.join(lines)
            base_name = job_def['jobName']
            if job_def['jobId'] in success_ids:
                base_name += '/SUCCESS'
            elif job_def['jobId'] in failure_ids:
                base_name += '/FAILED'
            else:
                logger.error("Job cannot be logged unless completed.")
                continue
            logger.info('Stashing ' + base_name)
            stash_log(log_str, base_name)
        except Exception as e:
            logger.error("Failed to save logs for: %s" % str(job_def_tpl))
            logger.exception(e)
    ids_stashed |= {jid for jids in [success_ids, failure_ids] for jid in jids}
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
                                "Got %d environments instead of 1."
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
    # Get the relevant instance ids from the ecs cluster
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
    # Get the Elsevier keys from the Elsevier client
    environment_vars = get_elsevier_api_keys()

    # Only include values that are not empty.
    return [var_dict for var_dict in environment_vars
            if var_dict['value'] and var_dict['name']]


class Submitter(object):
    _s3_input_name = NotImplemented
    _purpose = NotImplemented

    # The job queue on which these jobs will be submitted.
    _job_queue = NotImplemented

    # A dictionary of job_def names as keys, with a list of applicable readers
    # as the values.
    _job_def_dict = NotImplemented

    def __init__(self, basename, readers, project_name=None, **options):
        self.basename = basename
        if 'all' in readers:
            self.readers = [rc.name.lower() for rc in get_reader_classes()]
        else:
            self.readers = readers
        self.project_name = project_name
        self.job_list = None
        self.options = options
        self.ids_per_job = None
        self.running = None
        self.monitor = BatchMonitor(self._job_queue, self.job_list,
                                    self.basename)
        return

    def set_options(self, **kwargs):
        """Set the options of reading job."""
        # This should be more specifically implemented in a child class.
        self.options = kwargs
        return

    def _iter_commands(self, start_ix, end_ix):
        for job_def, reader_list in self._job_def_dict.items():
            reader_list = [r for r in reader_list if r in self.readers]
            if not reader_list:
                continue

            job_name = '%s_%d_%d' % (self.basename, start_ix, end_ix)
            job_name += '_' + '_'.join(reader_list)
            cmd = self._get_base(job_name, start_ix, end_ix)
            cmd += ['-r'] + reader_list
            cmd += self._get_extensions()
            for arg in cmd:
                if not isinstance(arg, str):
                    logger.warning("Argument of command is not a string: %s"
                               % repr(arg))
            yield job_name, cmd, job_def

    def _get_base(self, job_name, start_ix, end_ix):
        raise NotImplementedError

    def _get_extensions(self):
        return []

    def submit_reading(self, input_fname, start_ix, end_ix, ids_per_job,
                       num_tries=1, stagger=0):
        """Submit a batch of reading jobs

        Parameters
        ----------
        input_fname : str
            The name of the file containing the ids to be read.
        start_ix : int
            The line index of the first item in the list to read.
        end_ix : int
            The line index of the last item in the list to be read.
        ids_per_job : int
            The number of ids to be given to each job.
        num_tries : int
            The number of times a job may be attempted.
        stagger : float
            The number of seconds to wait between job submissions.

        Returns
        -------
        job_list : list[str]
            A list of job id strings.
        """
        self.job_list = []

        # stash this for later.
        self.ids_per_job = ids_per_job

        # Upload the pmid_list to Amazon S3
        id_list_key = 'reading_results/%s/%s' % (self.basename,
                                                 self._s3_input_name)
        s3_client = boto3.client('s3')
        s3_client.upload_file(input_fname, bucket_name, id_list_key)

        # If no end index is specified, read all the PMIDs
        if end_ix is None:
            with open(input_fname, 'rt') as f:
                lines = f.readlines()
                end_ix = len(lines)

        if start_ix is None:
            start_ix = 0

        # Get environment variables
        environment_vars = get_environment()

        # Iterate over the list of PMIDs and submit the job in chunks
        batch_client = boto3.client('batch', region_name='us-east-1')

        # Check to see if we've already been given a signal to quit.
        if self.running is None:
            self.running = True
        elif not self.running:
            return None

        for job_start_ix in range(start_ix, end_ix, ids_per_job):

            # Check for a stop signal
            if not self.running:
                logger.info("Running was switched off, discontinuing...")
                break

            # Generate the command for this batch.
            job_end_ix = job_start_ix + ids_per_job
            if job_end_ix > end_ix:
                job_end_ix = end_ix

            # Enter a command for reach job
            command_iter = self._iter_commands(job_start_ix, job_end_ix)
            for job_name, cmd, job_def in command_iter:
                command_list = get_batch_command(cmd, purpose=self._purpose,
                                                 project=self.project_name)
                logger.info('Command list: %s' % str(command_list))

                # Submit the job.
                job_info = batch_client.submit_job(
                    jobName=job_name,
                    jobQueue=self._job_queue,
                    jobDefinition=job_def,
                    containerOverrides={
                        'environment': environment_vars,
                        'command': command_list},
                    retryStrategy={'attempts': num_tries}
                )

                # Record the job id.
                logger.info("submitted...")
                self.job_list.append({'jobId': job_info['jobId']})
                logger.info("Sleeping for %d seconds..." % stagger)
                sleep(stagger)

        return self.job_list

    def watch_and_wait(self, poll_interval=10, idle_log_timeout=None,
                       kill_on_timeout=False, stash_log_method=None,
                       tag_instances=False, kill_on_exception=True, **kwargs):
        """This provides shortcut access to the wait_for_complete_function."""
        try:
            res = self.monitor.watch_and_wait(
                poll_interval=poll_interval, idle_log_timeout=idle_log_timeout,
                kill_on_log_timeout=kill_on_timeout,
                stash_log_method=stash_log_method, tag_instances=tag_instances,
                **kwargs
            )
        except (BaseException, KeyboardInterrupt) as e:
            logger.error("Exception in wait_for_complete:")
            logger.exception(e)
            if kill_on_exception:
                logger.info("Killing all my jobs...")
                kill_all(self._job_queue, kill_list=self.job_list,
                         reason='Exception in monitor, jobs aborted.')
            raise e
        return res

    def run(self, input_fname, ids_per_job, stagger=0, **wait_params):
        """Run this submission all the way.

        This method will run both `submit_reading` and `watch_and_wait`,
        blocking on the latter.
        """

        submit_thread = Thread(target=self.submit_reading,
                               args=(input_fname, 0, None, ids_per_job),
                               kwargs={'stagger': stagger},
                               daemon=True)
        submit_thread.start()
        try:
            logger.info("Waiting for just a sec...")
            sleep(1)
            wait_params['wait_for_first_job'] = True
            wait_params['kill_on_exception'] = True
            self.watch_and_wait(**wait_params)
            submit_thread.join(0)
            if submit_thread.is_alive():
                logger.warning("Submit thread is still running even after job "
                               "completion.")
        except BaseException as e:
            logger.error("Watch and wait failed...")
            logger.exception(e)
        finally:
            logger.info("Aborting jobs...")
            # Send a signal to the submission loop (on a thread) to stop.
            self.running = False
            submit_thread.join()
            print(submit_thread.is_alive())

        self.running = None
        return submit_thread


class PmidSubmitter(Submitter):
    _s3_input_name = 'pmids'
    _purpose = 'pmid_reading'
    _job_queue = 'run_reach_queue'
    _job_def_dict = {'run_reach_jobdef': ['reach', 'sparser'],
                     'run_db_reading_isi_jobdef': ['isi']}

    def _get_base(self, job_name, start_ix, end_ix):
        base = ['python', '-m', 'indra.tools.reading.pmid_reading.read_pmids_aws',
                self.basename, '/tmp', '16', str(start_ix), str(end_ix)]
        return base

    def _get_extensions(self):
        extensions = []
        for opt_key in ['force_read', 'force_fulltext']:
            if self.options.get(opt_key, False):
                extensions.append('--' + opt_key)
        return extensions

    def set_options(self, force_read=False, force_fulltext=False):
        """Set the options for this run."""
        self.options['force_read'] = force_read
        self.options['force_fulltext'] = force_fulltext
        return

    def submit_combine(self):
        job_ids = self.job_list
        if job_ids is not None and len(job_ids) > 20:
            print("WARNING: boto3 cannot support waiting for more than 20 jobs.")
            print("Please wait for the reading to finish, then run again with the")
            print("`combine` option.")
            return

        # Get environment variables
        environment_vars = get_environment()

        job_name = '%s_combine_reading_results' % self.basename
        command_list = get_batch_command(
            ['python', '-m', 'indra.tools.reading.assemble_reading_stmts_aws',
             self.basename, '-r'] + self.readers,
            purpose='pmid_reading',
            project=self.project_name
            )
        logger.info('Command list: %s' % str(command_list))
        kwargs = {'jobName': job_name, 'jobQueue': self._job_queue,
                  'jobDefinition': 'run_reach_jobdef',
                  'containerOverrides': {'environment': environment_vars,
                                         'command': command_list,
                                         'memory': 60000, 'vcpus': 1}}
        if job_ids:
            kwargs['dependsOn'] = job_ids
        batch_client = boto3.client('batch')
        batch_client.submit_job(**kwargs)
        logger.info("submitted...")
        return


def submit_reading(basename, pmid_list_filename, readers, start_ix=None,
                   end_ix=None, pmids_per_job=3000, num_tries=2,
                   force_read=False, force_fulltext=False, project_name=None):
    """Submit an old-style pmid-centered no-database s3 only reading job.

    This function is provided for the sake of backward compatibility. It is
    preferred that you use the object-oriented PmidSubmitter and the
    submit_reading job going forward.
    """
    sub = PmidSubmitter(basename, readers, project_name)
    sub.set_options(force_read, force_fulltext)
    sub.submit_reading(pmid_list_filename, start_ix, end_ix, pmids_per_job,
                       num_tries)
    return sub.job_list


def submit_combine(basename, readers, job_ids=None, project_name=None):
    """Submit a batch job to combine the outputs of a reading job.

    This function is provided for backwards compatibility. You should use the
    PmidSubmitter and submit_combine methods.
    """
    sub = PmidSubmitter(basename, readers, project_name)
    sub.job_list = job_ids
    sub.submit_combine()
    return sub


def create_submit_parser():
    import argparse
    parent_submit_parser = argparse.ArgumentParser(add_help=False)
    parent_submit_parser.add_argument(
        'basename',
        help='Defines job names and S3 keys'
    )
    parent_submit_parser.add_argument(
        '-r', '--readers',
        dest='readers',
        choices=[rc.name.lower() for rc in get_reader_classes()] + ['all'],
        default=['all'],
        nargs='+',
        help='Choose which reader(s) to use.'
    )
    parent_submit_parser.add_argument(
        '--project',
        help=('Set the project name. Default is DEFAULT_AWS_PROJECT in the '
              'config.')
    )
    return parent_submit_parser


def create_read_parser():
    import argparse
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
    return parent_read_parser


def create_parser():
    parent_submit_parser = create_submit_parser()
    parent_read_parser = create_read_parser()

    # Make non_db_parser and get subparsers
    parser = argparse.ArgumentParser(
        'indra.tools.reading.submit_reading_pipeline.py',
        description=('Run reading by collecting content, and save as pickles. '
                     'This option requires that ids are given as a list of '
                     'pmids, one line per pmid.'),
        epilog=('Note that `python wait_for_complete.py ...` should be run as '
                'soon as this command completes successfully. For more '
                'details use `python wait_for_complete.py -h`.')
    )
    subparsers = parser.add_subparsers(
        title='Job Type',
        help='Type of jobs to submit.'
    )
    subparsers.required = True
    subparsers.dest = 'job_type'

    # Create subparsers for the no-db option.
    subparsers.add_parser(
        'read',
        parents=[parent_read_parser, parent_submit_parser],
        help='Run reading on AWS batch and cache INDRA Statements on S3.',
        description='Run reading on batch and cache INDRA Statements on S3.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers.add_parser(
        'combine',
        parents=[parent_submit_parser],
        help='Combine INDRA Statement subsets into a single file.',
        description='Combine INDRA Statement subsets into a single file.'
    )
    subparsers.add_parser(
        'full',
        parents=[parent_read_parser, parent_submit_parser],
        help='Run reading and combine INDRA Statements when done.',
        description='Run reading and combine INDRA Statements when done.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    return parser


if __name__ == '__main__':

    parser = create_parser()
    args = parser.parse_args()

    job_ids = None
    sub = PmidSubmitter(args.basename, args.readers, args.project)
    sub.set_options(args.force_read, args.force_fulltext)
    if args.job_type in ['read', 'full']:
        sub.submit_reading(args.input_file, args.start_ix, args.end_ix,
                           args.ids_per_job)
    if args.job_type in ['combine', 'full']:
        sub.submit_combine()
