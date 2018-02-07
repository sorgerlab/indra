from argparse import ArgumentParser

if __name__ == '__main__':
    parser = ArgumentParser(
        'wait_for_complete.py',
        usage='%(prog)s [-h] queue_name [options]',
        description=('Wait for a set of batch jobs to complete, and monitor '
                     'them as they run.'),
        epilog=('Jobs can also be monitored, terminated, and otherwise '
                'managed on the AWS website. However this tool will also tag '
                'the instances, and should be run whenever a job is submitted '
                'to AWS.')
        )
    parser.add_argument(
        dest='queue_name',
        help=('The name of the queue to watch and wait for completion. If no '
              'jobs are specified, this will wait until all jobs in the queue '
              'are completed (either SUCCEEDED or FAILED).')
        )
    parser.add_argument(
        '--watch', '-w',
        dest='job_list',
        metavar='JOB_ID',
        nargs='+',
        help=('Specify particular jobs using their job ids, as reported by '
              'the submit command. Many ids may be specified.')
        )
    parser.add_argument(
        '--prefix', '-p',
        dest='job_name_prefix',
        help='Specify a prefix for the name of the jobs to watch and wait for.'
        )
    parser.add_argument(
        '--interval', '-i',
        dest='poll_interval',
        default=10,
        type=int,
        help=('The time interval to wait between job status checks, in '
              'seconds (default: %(default)d seconds).')
        )
    parser.add_argument(
        '--timeout', '-T',
        metavar='TIMEOUT',
        type=int,
        help=('If the logs are not updated for %(metavar)s seconds, '
              'print a warning. If `--kill_on_log_timeout` flag is set, then '
              'the offending jobs will be automatically terminated.')
        )
    parser.add_argument(
        '--kill_on_timeout', '-K',
        action='store_true',
        help='If a log times out, terminate the offending job.'
        )
    parser.add_argument(
        '--stash_log_method', '-l',
        choices=['s3', 'local'],
        metavar='METHOD',
        help=('Select a method from: [%(choices)s] to store the job logs. '
              'If no method is specified, the logs will not be '
              'loaded off of AWS. If \'s3\' is specified, then '
              '`job_name_prefix` must also be given, as this will indicate '
              'where on s3 to store the logs.')
        )
    args = parser.parse_args()

    from indra.tools.reading.submit_reading_pipeline import wait_for_complete

    job_list = None
    if args.job_list is not None:
        job_list = [{'jobId': jid} for jid in args.job_list]

    wait_for_complete(args.queue_name, job_list, args.job_name_prefix,
                      args.poll_interval, args.timeout,
                      args.kill_on_timeout, args.stash_log_method)
