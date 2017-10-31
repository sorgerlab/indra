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
    for status in ('STARTING', 'RUNNING'):
        running = batch.list_jobs(jobQueue=job_queue, jobStatus=status)
        job_info = running.get('jobSummaryList')
        if job_info:
            job_ids = [job['jobId'] for job in job_info]
            for job_id in job_ids:
                res = batch.terminate_job(jobId=job_id, reason=reason)


def get_jobs(job_queue='run_reach_queue', job_status='RUNNING'):
    """Returns a list of dicts with jobName and jobId for each job with the
    given status."""
    batch = boto3.client('batch')
    jobs = batch.list_jobs(jobQueue=job_queue, jobStatus=job_status)
    return jobs.get('jobSummaryList')


def get_job_log(job_info, log_group_name='/aws/batch/job',
                write_file=True):
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
    lines = []
    kwargs = {'logGroupName': log_group_name,
              'logStreamName': log_stream_name,
              'startFromHead': True}
    print("Getting log for %s/%s" % (job_name, job_id))
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
        print('%d %s' % (len(lines), lines[-1]))
    if write_file:
        with open('%s_%s.log' % (job_name, job_id), 'wt') as f:
            for line in lines:
                f.write(line)
    return lines


def dump_logs(job_queue='run_reach_queue', job_status='RUNNING'):
    """Write logs for all jobs with given the status to files."""
    jobs = get_jobs(job_queue, job_status)
    for job in jobs:
        get_job_log(job, write_file=True)

