
import re
import boto3
import pickle
import logging
from numpy import median, arange, array
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from datetime import datetime, timedelta

from indra.tools.reading.util.reporter import Reporter
from indra.util.get_version import get_git_info

from indra.util.nested_dict import NestedDict

from indra.tools.reading.submit_reading_pipeline import create_submit_parser, \
    create_read_parser, Submitter

bucket_name = 'bigmech'

logger = logging.getLogger('indra_db_reading')


class DbReadingSubmitter(Submitter):
    _s3_input_name = 'id_list'
    _purpose = 'db_reading'
    _job_queue = 'run_db_reading_queue'
    _job_def = 'run_db_reading_jobdef'

    def __init__(self, *args, **kwargs):
        super(DbReadingSubmitter, self).__init__(*args, **kwargs)
        self.time_tag = datetime.now().strftime('%Y%m%d_%H%M')
        self.reporter = Reporter(self.basename + '_summary_%s' % self.time_tag)
        self.reporter.sections = {'Plots': [], 'Totals': [], 'Git': []}
        self.reporter.set_section_order(['Git', 'Totals', 'Plots'])
        self.run_record = {}
        return

    def _get_base(self, job_name, start_ix, end_ix):
        read_mode = 'all' if self.options.get('force_read', False) else 'unread'
        stmt_mode = 'none' if self.options.get('no_stmts', False) else 'all'

        job_name = '%s_%d_%d' % (self.basename, start_ix, end_ix)
        base = ['python', '-m', 'indra.tools.reading.db_reading.read_db_aws',
                self.basename]
        base += [job_name]
        base += ['/tmp', read_mode, stmt_mode, '32', str(start_ix), str(end_ix)]
        return base

    def _get_extensions(self):
        extensions = []
        if self.options.get('force_fulltext', False):
            extensions.append('--force_fulltext')
        if self.options.get('prioritize', False):
            extensions.append('--read_best_fulltext')
        max_reach_input_len = self.options.get('max_reach_input_len')
        max_reach_space_ratio = self.options.get('max_reach_space_ratio')
        if max_reach_input_len is not None:
            extensions += ['--max_reach_input_len', max_reach_input_len]
        if max_reach_space_ratio is not None:
            extensions += ['--max_reach_space_ratio', max_reach_space_ratio]
        return extensions

    def set_options(self, force_read=False, no_stmts=False,
                    force_fulltext=False, prioritize=False,
                    max_reach_input_len=None, max_reach_space_ratio=None):
        self.options['force_fulltext'] = force_fulltext
        self.options['prioritize'] = prioritize
        self.options['max_reach_input_len'] = max_reach_input_len
        self.options['max_reach_space_ratio'] = max_reach_space_ratio
        return

    def watch_and_wait(self, *args, **kwargs):
        kwargs['result_record'] = self.run_record
        super(DbReadingSubmitter, self).watch_and_wait(*args, **kwargs)
        self.produce_report()

    @staticmethod
    def _parse_time(time_str):
        """Create a timedelta or datetime object from default string reprs."""
        try:
            # This is kinda terrible, but it is the easiest way to distinguish
            # them.
            if '-' in time_str:
                time_fmt = '%Y-%m-%d %H:%M:%S'
                if '.' in time_str:
                    pre_dec, post_dec = time_str.split('.')
                    dt = datetime.strptime(pre_dec, time_fmt)
                    dt.replace(microsecond=int(post_dec))
                else:
                    dt = datetime.strftime(time_str, time_fmt)
                return dt
            else:
                if 'day' in time_str:
                    m = re.match(('(?P<days>[-\d]+) day[s]*, '
                                  '(?P<hours>\d+):(?P<minutes>\d+):'
                                  '(?P<seconds>\d[\.\d+]*)'),
                                 time_str)
                else:
                    m = re.match(('(?P<hours>\d+):(?P<minutes>\d+):'
                                  '(?P<seconds>\d[\.\d+]*)'),
                                 time_str)
                return timedelta(**{key: float(val)
                                    for key, val in m.groupdict().items()})
        except Exception as e:
            logger.error('Failed to parse \"%s\".' % time_str)
            raise e

    def _get_results_file_tree(self, s3, s3_prefix):
        relevant_files = s3.list_objects(Bucket=bucket_name, Prefix=s3_prefix)
        file_tree = NestedDict()
        file_keys = [entry['Key'] for entry in relevant_files['Contents']]
        pref_path = s3_prefix.split('/')[:-1]   # avoid the trailing empty str.
        for key in file_keys:
            full_path = key.split('/')
            relevant_path = full_path[len(pref_path):]
            curr = file_tree
            for step in relevant_path:
                curr = curr[step]
            curr['key'] = key
        return file_tree

    def _get_txt_file_dict(self, file_bytes):
        line_list = file_bytes.decode('utf-8').splitlines()
        sc = ': '
        file_info = {}
        for line in line_list:
            segments = line.split(sc)
            file_info[segments[0].strip()] = sc.join(segments[1:]).strip()
        return file_info

    def _handle_git_info(self, ref, git_info, file_bytes):
        this_info = self._get_txt_file_dict(file_bytes)
        if git_info and this_info != git_info:
            logger.warning("Disagreement in git info in %s: "
                           "%s vs. %s."
                           % (ref, git_info, this_info))
        elif not git_info:
            git_info.update(this_info)
        return

    def _report_git_info(self, batch_git_info):
        self.reporter.add_text('Batch Git Info', section='Git', style='h1')
        for key, val in batch_git_info.items():
            label = key.replace('_', ' ').capitalize()
            self.reporter.add_text('%s: %s' % (label, val), section='Git')
        self.reporter.add_text('Launching System\'s Git Info', section='Git',
                               style='h1')
        git_info_dict = get_git_info()
        for key, val in git_info_dict.items():
            label = key.replace('_', ' ').capitalize()
            self.reporter.add_text('%s: %s' % (label, val), section='Git')
        return

    def _handle_timing(self, ref, timing_info, file_bytes):
        this_info = self._get_txt_file_dict(file_bytes)
        for stage, data in this_info.items():
            if stage not in timing_info.keys():
                logger.info("Adding timing stage: %s" % stage)
                timing_info[stage] = {}
            stage_info = timing_info[stage]
            timing_pairs = re.findall(r'(\w+):\s+([ 0-9:.\-]+)', data)
            if len(timing_pairs) is not 3:
                logger.warning("Not all timings present for %s "
                               "in %s." % (stage, ref))
            for label, time_str in timing_pairs:
                if label not in stage_info.keys():
                    stage_info[label] = {}
                # e.g. timing_info['reading']['start']['job_name'] = <datetime>
                stage_info[label][ref] = self._parse_time(time_str)
        return

    def _report_timing(self, timing_info):
        # Pivot the timing info.
        idx_patt = re.compile('%s_(\d+)_(\d+)' % self.basename)
        job_segs = NestedDict()
        plot_set = set()
        for stage, stage_d in timing_info.items():
            # e.g. reading, statement production...
            for metric, metric_d in stage_d.items():
                # e.g. start, end, ...
                for job_name, t in metric_d.items():
                    # e.g. job_basename_startIx_endIx
                    job_segs[job_name][stage][metric] = t
                    m = idx_patt.match(job_name)
                    if m is None:
                        logger.error("Unexpectedly formatted name: %s."
                                     % job_name)
                        continue
                    key = tuple([int(n) for n in m.groups()] + [job_name])
                    plot_set.add(key)
        plot_list = list(plot_set)
        plot_list.sort()

        # Use this for getting the minimum and maximum.
        all_times = [dt for job in job_segs.values() for stage in job.values()
                     for metric, dt in stage.items() if metric != 'duration']
        all_start = min(all_times)
        all_end = max(all_times)

        def get_time_tuple(stage_data):
            start_seconds = (stage_data['start'] - all_start).total_seconds()
            return start_seconds, stage_data['duration'].total_seconds()

        # Make the broken barh plots.
        w = 6.5
        h = 9
        fig = plt.figure(figsize=(w, h))
        gs = plt.GridSpec(2, 1, height_ratios=[10, 1])
        ax0 = plt.subplot(gs[0])
        ytick_pairs = []
        stages = ['reading', 'statement production', 'stats']
        t = arange((all_end - all_start).total_seconds())
        counts = dict.fromkeys(['jobs'] + stages)
        for k in counts.keys():
            counts[k] = array([0 for _ in t])
        for i, job_tpl in enumerate(plot_list):
            s_ix, e_ix, job_name = job_tpl
            job_d = job_segs[job_name]
            xs = [get_time_tuple(job_d[stg]) for stg in stages]
            ys = (s_ix, (e_ix - s_ix)*0.9)
            ytick_pairs.append(((s_ix + e_ix)/2, '%s_%s' % (s_ix, e_ix)))
            logger.debug("Making plot for: %s" % str((job_name, xs, ys)))
            ax0.broken_barh(xs, ys, facecolors=('red', 'green', 'blue'))

            for n, stg in enumerate(stages):
                cs = counts[stg]
                start = xs[n][0]
                dur = xs[n][1]
                cs[(t>start) & (t<(start + dur))] += 1
            cs = counts['jobs']
            cs[(t>xs[0][0]) & (t<(xs[-1][0] + xs[-1][1]))] += 1

        # Format the plot
        ax0.tick_params(top='off', left='off', right='off', bottom='off',
                       labelleft='on', labelbottom='off')
        for spine in ax0.spines.values():
            spine.set_visible(False)
        total_time = (all_end - all_start).total_seconds()
        ax0.set_xlim(0, total_time)
        ax0.set_ylabel(self.basename + '_ ...')
        print(ytick_pairs)
        yticks, ylabels = zip(*ytick_pairs)
        print(yticks)
        if not self.ids_per_job:
            print([yticks[i+1] - yticks[i]
                   for i in range(len(yticks) - 1)])
            # Infer if we don't have it.
            spacing = median([yticks[i+1] - yticks[i]
                              for i in range(len(yticks) - 1)])
            spacing = max(1, spacing)
        else:
            spacing = self.ids_per_job
        print(spacing)
        print(yticks[0], yticks[-1])
        ytick_range = list(arange(yticks[0], yticks[-1] + spacing, spacing))
        ylabel_filled = []
        for ytick in ytick_range:
            if ytick in yticks:
                ylabel_filled.append(ylabels[yticks.index(ytick)])
            else:
                ylabel_filled.append('FAILED')
        ax0.set_ylim(0, max(ytick_range) + spacing)
        ax0.set_yticks(ytick_range)
        ax0.set_yticklabels(ylabel_filled)

        # Plot the lower axis.
        legend_list = []
        color_map = {'jobs': 'k', 'reading': 'r', 'statement production': 'g',
                     'stats': 'b'}
        ax1 = plt.subplot(gs[1], sharex=ax0)
        for k, cs in counts.items():
            legend_list.append(k)
            ax1.plot(t, cs, color=color_map[k])
        for lbl, spine in ax1.spines.items():
            spine.set_visible(False)
        max_n = max(counts['jobs'])
        ax1.set_ylim(0, max_n + 1)
        ax1.set_xlim(0, total_time)
        yticks = list(range(0, max_n-max_n//5, max(1, max_n//5)))
        ax1.set_yticks(yticks + [max_n])
        ax1.set_yticklabels([str(n) for n in yticks] + ['max=%d' % max_n])
        ax1.set_ylabel('N_jobs')
        ax1.set_xlabel('Time since beginning [seconds]')

        # Make the figure borders more sensible.
        fig.tight_layout()
        img_path = 'time_figure.png'
        fig.savefig(img_path)
        self.reporter.add_image(img_path, width=w, height=h, section='Plots')
        return

    def _handle_sum_data(self, job_ref, summary_info, file_bytes):
        one_sum_data_dict = pickle.loads(file_bytes)
        for k, v in one_sum_data_dict.items():
            if k not in summary_info.keys():
                summary_info[k] = {}
            summary_info[k][job_ref] = v
        return

    def _report_sum_data(self, summary_info):
        # Two kind of things to handle:
        for k, job_dict in summary_info.items():
            if isinstance(list(job_dict.values())[0], dict):
                continue

            # Overall totals
            self.reporter.add_text('total %s: %d' % (k, sum(job_dict.values())),
                                   section='Totals')

            # Hists of totals.
            if len(job_dict) <= 1:
                continue

            w = 6.5
            h = 4
            fig = plt.figure(figsize=(w, h))
            plt.hist(list(job_dict.values()), align='left')
            plt.xlabel(k)
            plt.ylabel('Number of Jobs')
            fig.tight_layout()
            fname = k + '_hist.png'
            fig.savefig(fname)
            self.reporter.add_image(fname, width=w, height=h, section='Plots')
        return

    def _handle_hist_data(self, job_ref, hist_dict, file_bytes):
        a_hist_data_dict = pickle.loads(file_bytes)
        for k, v in a_hist_data_dict.items():
            if k not in hist_dict.keys():
                hist_dict[k] = {}
            hist_dict[k][job_ref] = v
        return

    def _report_hist_data(self, hist_dict):
        for k, data_dict in hist_dict.items():
            w = 6.5
            if k == ('stmts', 'readers'):
                h = 6
                fig = plt.figure(figsize=(w, h))
                data = {}
                for job_datum in data_dict.values():
                    for rdr, num in job_datum['data'].items():
                        if rdr not in data.keys():
                            data[rdr] = [num]
                        else:
                            data[rdr].append(num)
                N = len(data)
                key_list = list(data.keys())
                xtick_locs = arange(N)
                n = (N+1)*100 + 11
                ax0 = plt.subplot(n)
                ax0.bar(xtick_locs, [sum(data[k]) for k in key_list],
                        align='center')
                ax0.set_xticks(xtick_locs, key_list)
                ax0.set_xlabel('readers')
                ax0.set_ylabel('stmts')
                ax0.set_title('Reader production')
                rdr_ax_list = []
                for rdr, stmt_counts in data.items():
                    n += 1
                    if not rdr_ax_list:
                        ax = plt.subplot(n)
                    else:
                        ax = plt.subplot(n, sharex=rdr_ax_list[0])
                    ax.set_title(rdr)
                    ax.hist(stmt_counts, align='left')
                    ax.set_ylabel('jobs')
                    rdr_ax_list.append(ax)
                if rdr_ax_list:
                    ax.set_xlabel('stmts')
            else:  # TODO: Handle other summary plots.
                continue
            figname = '_'.join(k) + '.png'
            fig.savefig(figname)
            self.reporter.add_image(figname, width=w, height=h, section='Plots')

        return

    def produce_report(self):
        """Produce a report of the batch jobs."""
        s3_prefix = 'reading_results/%s/logs/%s/' % (self.basename,
                                                     self._job_queue)
        logger.info("Producing batch report for %s, from prefix %s."
                    % (self.basename, s3_prefix))
        s3 = boto3.client('s3')
        file_tree = self._get_results_file_tree(s3, s3_prefix)
        logger.info("Found %d relevant files." % len(file_tree))
        stat_files = {
            'git_info.txt': (self._handle_git_info, self._report_git_info),
            'timing.txt': (self._handle_timing, self._report_timing),
            'raw_tuples.pkl': (None, None),
            'hist_data.pkl': (self._handle_hist_data, self._report_hist_data),
            'sum_data.pkl': (self._handle_sum_data, self._report_sum_data)
            }
        stat_aggs = {}
        for stat_file, (handle_stats, report_stats) in stat_files.items():
            logger.info("Aggregating %s..." % stat_file)
            # Prep the data storage.
            my_agg = {}

            # Get a list of the relevant files (one per job).
            file_paths = file_tree.get_paths(stat_file)
            logger.info("Found %d files for %s." % (len(file_paths), stat_file))

            # Aggregate the data from all the jobs for each file type.
            for sub_path, file_entry in file_paths:
                s3_key = file_entry['key']
                ref = sub_path[0]
                file = s3.get_object(Bucket=bucket_name, Key=s3_key)
                file_bytes = file['Body'].read()
                if handle_stats is not None:
                    handle_stats(ref, my_agg, file_bytes)

            if report_stats is not None and len(my_agg):
                report_stats(my_agg)

            stat_aggs[stat_file] = my_agg

        for end_type, jobs in self.run_record.items():
            self.reporter.add_text('Jobs %s: %d' % (end_type, len(jobs)),
                                   section='Totals')

        s3_prefix = 'reading_results/%s/' % self.basename
        fname = self.reporter.make_report()
        with open(fname, 'rb') as f:
            s3.put_object(Bucket=bucket_name,
                          Key= s3_prefix + fname,
                          Body=f.read())
        s3.put_object(Bucket=bucket_name,
                      Key=s3_prefix + 'stat_aggregates_%s.pkl' % self.time_tag,
                      Body=pickle.dumps(stat_aggs))
        return file_tree, stat_aggs


def submit_db_reading(basename, id_list_filename, readers, start_ix=None,
                      end_ix=None, pmids_per_job=3000, num_tries=2,
                      force_read=False, force_fulltext=False,
                      read_all_fulltext=False, project_name=None,
                      max_reach_input_len=None, max_reach_space_ratio=None,
                      no_stmts=False):
    """Submit batch reading jobs that uses the database for content and results.

    This function is provided for backwards compatibility, use DbReadingSubmitter
    and its submit_reading method instead.
    """
    sub = DbReadingSubmitter(basename, readers, project_name)
    sub.set_options(force_read, no_stmts, force_fulltext, read_all_fulltext,
                    max_reach_input_len, max_reach_space_ratio)
    sub.submit_reading(id_list_filename, start_ix, end_ix, pmids_per_job,
                       num_tries)
    return sub


if __name__ == '__main__':
    import argparse

    parent_submit_parser = create_submit_parser()
    parent_read_parser = create_read_parser()

    # Make db parser and get subparsers.
    parser = argparse.ArgumentParser(
        'indra_db.reading.submig_reading_pipeline.py',
        parents=[parent_submit_parser, parent_read_parser],
        description=('Run reading with content on the db and submit results. '
                     'In this option, ids in \'input_file\' are given in the '
                     'format \'<id type>:<id>\'. Unlike no-db, there is no '
                     'need to combine pickles, and therefore no need to '
                     'specify your task further.'),
        )
    '''Not currently supported
    parent_db_parser.add_argument(
        '--no_upload',
        action='store_true',
        help='Don\'t upload results to the database.'
        )
    '''
    parser.add_argument(
        '--read_best_fulltext',
        action='store_true',
        help='Read only the best fulltext for input ids.'
    )
    parser.add_argument(
        '--no_statements',
        action='store_true',
        help='Choose to not produce any Statements; only readings will be done.'
    )
    parser.add_argument(
        '--max_reach_space_ratio',
        type=float,
        help='Set the maximum ratio of spaces to non-spaces for REACH input.',
        default=None
    )
    parser.add_argument(
        '--max_reach_input_len',
        type=int,
        help='Set the maximum length of content that REACH will read.',
        default=None
    )
    args = parser.parse_args()

    sub = DbReadingSubmitter(args.basename, args.readers, args.project)
    sub.set_options(args.force_read, args.no_statements,
                    args.force_fulltext, args.prioritize,
                    args.max_reach_input_len, args.max_reach_space_ratio)
    sub.submit_reading(args.input_file, args.start_ix, args.end_ix,
                       args.ids_per_job)
