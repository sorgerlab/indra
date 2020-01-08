def get_s3_root(s3_base, job_queue="run_db_reading_queue"):
    return s3_base + 'logs/%s/' % job_queue


def get_s3_job_prefix(s3_base, job_name, job_queue="run_db_reading_queue"):
    s3_root = get_s3_root(s3_base, job_queue)
    return s3_root + '%s/' % job_name


def get_s3_and_job_prefixes(base_name, group_name=None):
    """Returns the s3 prefix and job prefix."""
    if not group_name:
        s3_base = base_name
        job_base = base_name
    else:
        s3_base, job_base = [group_name + d + base_name for d in ['/', '_']]
    return 'reading_results/' + s3_base + '/', job_base
