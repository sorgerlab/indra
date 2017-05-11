import boto3

def kill_all(job_queue, reason='None given'):
    import ipdb; ipdb.set_trace()
    batch = boto3.client('batch')
    runnable = batch.list_jobs(jobQueue=job_queue, jobStatus='RUNNABLE')
    job_info = runnable.get('jobSummaryList')
    if job_info:
        job_ids = [job['jobId'] for job in job_info]
        # Cancel jobs
        for job_id in job_ids:
            batch.cancel_job(jobId=job_id, reason=reason)
    running = batch.list_jobs(jobQueue=job_queue, jobStatus='RUNNING')
    job_info = running.get('jobSummaryList')
    if job_info:
        job_ids = [job['jobId'] for job in job_info]
        for job_id in job_ids:
            batch.terminate_job(jobId=job_id, reason=reason)

