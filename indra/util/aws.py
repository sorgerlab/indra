import re
import boto3

def kill_all(job_queue, reason='None given'):
    """Terminates/cancels all RUNNING, RUNNABLE, and STARTING jobs."""
    batch = boto3.client('batch')
    runnable = batch.list_jobs(jobQueue=job_queue, jobStatus='RUNNABLE')
    job_info = runnable.get('jobSummaryList')
    if job_info:
        job_ids = [job['jobId'] for job in job_info]
        # Cancel jobs
        for job_id in job_ids:
            batch.cancel_job(jobId=job_id, reason=reason)
    for status in ('STARTING', 'RUNNABLE', 'RUNNING'):
        running = batch.list_jobs(jobQueue=job_queue, jobStatus=status)
        job_info = running.get('jobSummaryList')
        if job_info:
            job_ids = [job['jobId'] for job in job_info]
            for job_id in job_ids:
                print('Killing %s' % job_id)
                res = batch.terminate_job(jobId=job_id, reason=reason)


def get_jobs(job_queue='run_reach_queue', job_status='RUNNING'):
    """Returns a list of dicts with jobName and jobId for each job with the
    given status."""
    batch = boto3.client('batch')
    jobs = batch.list_jobs(jobQueue=job_queue, jobStatus=job_status)
    return jobs.get('jobSummaryList')


def get_job_log(job_info, log_group_name='/aws/batch/job',
                write_file=True, verbose=False):
    """Gets the Cloudwatch log associated with the given job.

    Parameters
    ----------
    job_info : dict
        dict containing entries for 'jobName' and 'jobId', e.g., as returned
        by get_jobs()
    log_group_name : string
        Name of the log group; defaults to '/aws/batch/job'
    write_file : boolean
        If True, writes the downloaded log to a text file with the filename
        '%s_%s.log' % (job_name, job_id)


    Returns
    -------
    list of strings
        The event messages in the log, with the earliest events listed first.
    """
    job_name = job_info['jobName']
    job_id = job_info['jobId']
    logs = boto3.client('logs')
    batch = boto3.client('batch')
    job_description = batch.describe_jobs(jobs=[job_id])
    log_stream_name = job_description['jobs'][0]['container']['logStreamName']
    stream_resp = logs.describe_log_streams(
                            logGroupName=log_group_name,
                            logStreamNamePrefix=log_stream_name)
    streams = stream_resp.get('logStreams')
    if not streams:
        print('No streams for job')
        return None
    elif len(streams) > 1:
        print('More than 1 stream for job, returning first')
    log_stream_name = streams[0]['logStreamName']
    if verbose:
        print("Getting log for %s/%s" % (job_name, job_id))
    out_file = ('%s_%s.log' % (job_name, job_id)) if write_file else None
    lines = get_log_by_name(log_group_name, log_stream_name, out_file, verbose)
    return lines


def get_log_by_name(log_group_name, log_stream_name, out_file=None,
                    verbose=True):
    """Download a log given the log's group and stream name.

    Parameters
    ----------
    log_group_name : str
        The name of the log group, e.g. /aws/batch/job.

    log_stream_name : str
        The name of the log stream, e.g. run_reach_jobdef/default/<UUID>

    Returns
    -------
    lines : list[str]
        The lines of the log as a list.
    """
    logs = boto3.client('logs')
    kwargs = {'logGroupName': log_group_name,
              'logStreamName': log_stream_name,
              'startFromHead': True}
    lines = []
    while True:
        response = logs.get_log_events(**kwargs)
        # If we've gotten all the events already, the nextForwardToken for
        # this call will be the same as the last one
        if response.get('nextForwardToken') == kwargs.get('nextToken'):
            break
        else:
            events = response.get('events')
            if events:
                lines += ['%s: %s\n' % (evt['timestamp'], evt['message'])
                          for evt in events]
            kwargs['nextToken'] = response.get('nextForwardToken')
        if verbose:
            print('%d %s' % (len(lines), lines[-1]))
    if out_file:
        with open(out_file, 'wt') as f:
            for line in lines:
                f.write(line)
    return lines


def dump_logs(job_queue='run_reach_queue', job_status='RUNNING'):
    """Write logs for all jobs with given the status to files."""
    jobs = get_jobs(job_queue, job_status)
    for job in jobs:
        get_job_log(job, write_file=True)


def analyze_reach_log(log):
    """Return unifinished PMIDs given a log file name."""
    def get_content_nums(txt):
        pat = 'Retrieved content for ([\d]+) / ([\d]+) papers to be read'
        res = re.match(pat, txt)
        has_content, total = res.groups() if res else None, None
        return has_content, total

    def get_started(txt):
        pat = 'Starting ([\d]+)'
        pmids = re.findall(pat, txt)
        return pmids

    def get_finished(txt):
        # TODO: it might be interesting to get the time it took to read
        # each paper here
        pat = 'Finished ([\d]+)'
        pmids = re.findall(pat, txt)
        return pmids

    with open(log, 'r') as fh:
        log = fh.read()
    has_content, total = get_content_nums(log)
    pmids_started = get_started(log)
    pmids_finished = get_finished(log)
    pmids_not_done = set(pmids_started) - set(pmids_finished)
    return pmids_not_done

