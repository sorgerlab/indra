import re
import boto3
from os.path import join


def analyze_reach_log(log_fname=None, log_str=None):
    """Return unifinished PMIDs given a log file name."""
    assert bool(log_fname) ^ bool(log_str), 'Must specify log_fname OR log_str'
    started_patt = re.compile('Starting ([\d]+)')
    # TODO: it might be interesting to get the time it took to read
    # each paper here
    finished_patt = re.compile('Finished ([\d]+)')

    def get_content_nums(txt):
        pat = 'Retrieved content for ([\d]+) / ([\d]+) papers to be read'
        res = re.match(pat, txt)
        has_content, total = res.groups() if res else None, None
        return has_content, total

    if log_fname:
        with open(log_fname, 'r') as fh:
            log_str = fh.read()
    # has_content, total = get_content_nums(log_str)  # unused
    pmids = {}
    pmids['started'] = started_patt.findall(log_str)
    pmids['finished'] = finished_patt.findall(log_str)
    pmids['not_done'] = set(pmids['started']) - set(pmids['finished'])
    return pmids


#==============================================================================
# Functions for analyzing a db reading submission
#==============================================================================


def get_logs_from_db_reading(job_prefix, reading_queue='run_db_reading_queue'):
    """Get the logs stashed on s3 for a particular reading."""
    s3 = boto3.client('s3')
    gen_prefix = 'reading_results/%s/logs/%s' % (job_prefix, reading_queue)
    job_log_data = s3.list_objects_v2(Bucket='bigmech',
                                      Prefix=join(gen_prefix, job_prefix))
    # TODO: Track success/failure
    log_strs = []
    for fdict in job_log_data['Contents']:
        resp = s3.get_object(Bucket='bigmech', Key=fdict['Key'])
        log_strs.append(resp['Body'].read().decode('utf-8'))
    return log_strs


def separate_reach_logs(log_str):
    """Get the list of reach logs from the overall logs."""
    log_lines = log_str.splitlines()
    reach_logs = []
    reach_lines = []
    adding_reach_lines = False
    for l in log_lines[:]:
        if not adding_reach_lines and 'Beginning reach' in l:
            adding_reach_lines = True
        elif adding_reach_lines and 'Reach finished' in l:
            adding_reach_lines = False
            reach_logs.append(('SUCCEEDED', '\n'.join(reach_lines)))
            reach_lines = []
        elif adding_reach_lines:
            reach_lines.append(l.split('readers - ')[1])
            log_lines.remove(l)
    if adding_reach_lines:
        reach_logs.append(('FAILURE', '\n'.join(reach_lines)))
    return '\n'.join(log_lines), reach_logs


def get_top_level_summary_of_log(log_str):
    ret_str = 'Event Summary:'
    ret_str += '\n' + '-'*len(ret_str)
    ret_str += '\nUseful INFO:\n  '
    ret_str += '\n  '.join(re.findall(
        ('INFO: \[.*?\] indra/((?!readers).* - '
         '(?!Got no statements|Saving sparser)(?=.*\d.*).*)'),
        log_str))
    ret_str += '\nWARNINGS that occured:\n  '
    ret_str += '\n  '.join(set(get_indra_logs_by_priority(log_str, 'WARNING')))
    ret_str += '\nERRORS that occured:\n  '
    ret_str += '\n  '.join(set(get_indra_logs_by_priority(log_str, 'ERROR')))
    return ret_str


def get_top_level_summary_of_logs(log_str_list):
    ret_dict = {}
    ret_dict['total_stats'] = {}
    ret_dict['err_set'] = set()
    ret_dict['warn_set'] = set()
    ret_dict['unyielding_tcids'] = set()
    ret_dict['num_failures'] = 0
    for log_str in log_str_list:
        try:
            stat_dict = get_reading_stats(log_str)
            ret_dict['total_stats'] = {k: ret_dict['total_stats'].get(k, 0) + v
                                       for k, v in stat_dict.items()}
        except GetReadingStatsError:
            ret_dict['num_failures'] += 1
        ret_dict['err_set'] |= set(get_indra_logs_by_priority(log_str,
                                                              'ERROR'))
        ret_dict['warn_set'] |= set(get_indra_logs_by_priority(log_str,
                                                               'WARNING'))
        ret_dict['unyielding_tcids'] |= get_unyielding_tcids(log_str)
    ret_dict['err_tcids'] = {int(re.findall('(\d+)', err_str)[0])
                             for err_str in ret_dict['err_set']
                             if 'Got exception creating statements' in err_str}
    return ret_dict


def get_indra_logs_by_priority(log_str, priority='INFO'):
    return re.findall('%s: \[.*?\] indra/(.*)' % priority, log_str)


def get_unyielding_tcids(log_str):
    """Extract the set of tcids for which no statements were created."""
    tcid_strs = re.findall('INFO: \[.*?\].*? - Got no statements for (\d+).*',
                           log_str)
    return {int(tcid_str) for tcid_str in tcid_strs}


class GetReadingStatsError(Exception):
    pass


def get_reading_stats(log_str):
    def re_get_nums(patt_str, default=None):
        re_ret = re.search(patt_str, log_str)
        if re_ret is not None:
            nums = [int(num_str) for num_str in re_ret.groups()]
        elif default is None:
            raise GetReadingStatsError("couldn't match patt \"%s\"" % patt_str)
        else:
            nums = [default]*patt_str.count('(\d+)')
        return nums
    ret_dict = {}
    ret_dict['num_prex_readings'] = \
        re_get_nums('Found (\d+) pre-existing readings', 0)[0]
    try:
        ret_dict['num_new_readings'] = re_get_nums('Made (\d+) new readings')[0]
    except:
        ret_dict['num_new_readings'] = None
    ret_dict['num_succeeded'] = \
        re_get_nums('Adding (\d+)/\d+ reading entries')[0]
    ret_dict['num_stmts'], ret_dict['num_readings'] = \
        re_get_nums('Found (\d+) statements from (\d+) readings')
    ret_dict['num_agents'] = \
        re_get_nums('Received request to copy (\d+) entries '
                    'into .{3,4}agents')[0]
    ret_dict['num_statements'] = \
        re_get_nums('Received request to copy (\d+) entries into '
                    '.{3,4}statements')[0]
    return ret_dict


def analyze_db_reading(job_prefix, reading_queue='run_db_reading_queue'):
    """Run various analysis on a particular reading job."""
    # Analyze reach failures
    log_strs = get_logs_from_db_reading(job_prefix, reading_queue)
    indra_log_strs = []
    all_reach_logs = []
    log_stats = []
    for log_str in log_strs:
        log_str, reach_logs = separate_reach_logs(log_str)
        all_reach_logs.extend(reach_logs)
        indra_log_strs.append(log_str)
        log_stats.append(get_reading_stats(log_str))

    # Analayze the reach failures.
    failed_reach_logs = [reach_log_str
                         for result, reach_log_str in all_reach_logs
                         if result == 'FAILURE']
    failed_id_dicts = [analyze_reach_log(log_str=reach_log)
                       for reach_log in failed_reach_logs if bool(reach_log)]
    tcids_unfinished = {id_dict['not_done'] for id_dict in failed_id_dicts}
    print("Found %d unfinished tcids." % len(tcids_unfinished))

    # Summarize the global stats
    if log_stats:
        sum_dict = dict.fromkeys(log_stats[0].keys())
        for log_stat in log_stats:
            for k in log_stat.keys():
                if isinstance(log_stat[k], list):
                    if k not in sum_dict.keys():
                        sum_dict[k] = [0]*len(log_stat[k])
                    sum_dict[k] = [sum_dict[k][i] + log_stat[k][i]
                                   for i in range(len(log_stat[k]))]
                else:
                    if k not in sum_dict.keys():
                        sum_dict[k] = 0
                    sum_dict[k] += log_stat[k]
    else:
        sum_dict = {}

    return tcids_unfinished, sum_dict, log_stats
