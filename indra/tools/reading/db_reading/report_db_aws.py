from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import pickle
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from datetime import datetime

from indra.util.get_version import get_git_info
from indra.tools.reading.util.reporter import Reporter


class DbAwsStatReporter(Reporter):
    """A class to handle generating the reports made at the end of a job.

    Parameters
    ----------
    job_name : str
        The name of the job about which this report is being produced.
    s3_log_prefix: str
        The desired prefix location in which the information about the job is
        stored on s3.
    s3 : boto3.client('s3') instance
        An s3 client from boto3.
    bucket_name : str
        The name of the root s3 bucket.
    """
    def __init__(self, job_name, s3_log_prefix, s3, bucket_name):
        super(DbAwsStatReporter, self).__init__('%s_summary' % job_name)
        self.s3_prefix = s3_log_prefix + 'statistics/'
        self.bucket_name = bucket_name
        self.s3 = s3
        self.summary_dict = {}
        self.hist_dict = {}
        self.set_title("Report of Database Reading Batch Job")
        self.sections = {'Summary Statistics': [], 'Plots': [], 'Git Info': [],
                         'Job Info': []}
        self._make_job_line('Job Name', job_name)
        self._make_job_line('Job s3 prefix', s3_log_prefix)
        self._get_git_info()
        self.set_section_order(['Job Info', 'Git Info', 'Summary Statistics',
                                'Plots'])
        return

    def _get_git_info(self):
        """Stash the git info in a text file."""
        git_info_dict = get_git_info()
        text_file_content = ''
        for key, val in git_info_dict.items():
            label = key.replace('_', ' ').capitalize()
            text_file_content += '%s: %s\n' % (label, val)
            self.add_text('%s: %s' % (label, val), section='Git Info')
        self.s3.put_object(Key=self.s3_prefix + 'git_info.txt',
                           Body=text_file_content, Bucket=self.bucket_name)
        return

    def _plot_hist(self, agged, agg_over, data):
        """Make and stash a figure for histogram-like data."""
        fig = plt.figure()

        # Check if it's worth actually making a histogram...
        if not len(data):
            self.add_text('No data for %s vs %s.' % (agged, agg_over),
                          section='Plots')
            return
        if len(set(data)) == 1:
            if isinstance(data, dict):
                k = list(data.keys())[0]
                v = list(data.values())[0]
                txt = 'One of the %s, %s, has %d %s.' % (agg_over, k, v, agged)
                self.add_text(txt, section='Plots')
            else:
                self.add_text('%s per %s: %d' % (agged, agg_over, data[0]),
                              section='Plots')
            return

        # Make the histogram.
        if isinstance(data, dict):
            key_list = list(data.keys())
            xtick_locs = np.arange(len(data))
            plt.bar(xtick_locs, [data[k] for k in key_list],
                     align='center')
            plt.xticks(xtick_locs, key_list)
            plt.xlabel(agg_over)
            plt.ylabel(agged)
        else:
            plt.hist(data, bins=np.arange(min(data), max(data)), log=True,
                     align='left', edgecolor='none')
            plt.xlim(min(data), max(data))
            plt.xlabel(agged)
            plt.ylabel(agg_over)
        fname = '%s_per_%s.png' % (agged, agg_over)
        fig.set_size_inches(6, 4)
        fig.tight_layout()
        fig.savefig(fname)
        with open(fname, 'rb') as f:
            s3_key = self.s3_prefix + fname
            self.s3.put_object(Key=s3_key, Body=f.read(),
                               Bucket=self.bucket_name)
        self.add_image(fname, width=6, height=4, section='Plots')
        return

    def _make_timing_report(self, starts, ends):
        """Stash a text file with the timings: start, end, and duration."""
        # Report on the timing
        timing_str = ''
        for step in ['reading', 'statement production', 'stats']:
            time_taken = ends[step] - starts[step]
            timing_str += ('%22s: start: %s, end: %s, duration: %s\n'
                           % (step, starts[step], ends[step], time_taken))

        self.s3.put_object(Key=self.s3_prefix + 'timing.txt', Body=timing_str,
                           Bucket=self.bucket_name)
        return

    def _stash_data(self):
        """Store the data in pickle files. This should be done last."""
        self.s3.put_object(Key=self.s3_prefix + 'hist_data.pkl',
                           Body=pickle.dumps(self.hist_dict),
                           Bucket=self.bucket_name)
        self.s3.put_object(Key=self.s3_prefix + 'sum_data.pkl',
                           Bucket=self.bucket_name,
                           Body=pickle.dumps(self.summary_dict))
        return

    def _populate_hist_data(self, readings_with_stmts, readings_with_no_stmts):
        """Aggregate data from the various readings into arrays."""
        # Do a bunch of aggregation
        tc_rd_dict = {}
        tc_stmt_dict = {}
        rd_stmt_dict = {}
        reader_stmts = {}
        reader_tcids = {}
        reader_rids = {}
        for rid, tcid, reader, stmts in readings_with_stmts:
            # Handle things keyed by tcid
            if tcid not in tc_rd_dict.keys():
                tc_rd_dict[tcid] = {rid}
                tc_stmt_dict[tcid] = set(stmts)
            else:
                tc_rd_dict[tcid].add(rid)
                tc_stmt_dict[tcid] |= set(stmts)

            # Handle things keyed by rid
            if rid not in rd_stmt_dict.keys():
                rd_stmt_dict[rid] = set(stmts)
            else:
                rd_stmt_dict[rid] |= set(stmts)  # this shouldn't really happen.

            # Handle things keyed by reader
            if reader not in reader_stmts.keys():
                reader_stmts[reader] = set(stmts)
                reader_tcids[reader] = {tcid}
                reader_rids[reader] = {rid}
            else:
                reader_stmts[reader] |= set(stmts)
                reader_tcids[reader].add(tcid)
                reader_rids[reader].add(rid)

        for rid, tcid, reader, _ in readings_with_no_stmts:
            # Handle things keyed by tcid
            if tcid not in tc_rd_dict.keys():
                tc_rd_dict[tcid] = {rid}
            else:
                tc_rd_dict[tcid].add(rid)

            # Handle things keyed by reader
            if reader not in reader_tcids.keys():
                reader_tcids[reader] = {tcid}
                reader_rids[reader] = {rid}
            else:
                reader_tcids[reader].add(tcid)
                reader_rids[reader].add(rid)

        # Produce some numpy count arrays.
        self.hist_dict[('readings', 'text content')] = {
            'data': np.array([len(rid_set) for rid_set in tc_rd_dict.values()])
            }
        self.hist_dict[('stmts', 'text content')] = {
            'data': np.array([len(stmts) for stmts in tc_stmt_dict.values()])
            }
        self.hist_dict[('stmts', 'readings')] = {
            'data': np.array([len(stmts) for stmts in rd_stmt_dict.values()])
            }
        self.hist_dict[('stmts', 'readers')] = {
            'data': {rdr: len(stmts) for rdr, stmts in reader_stmts.items()}
            }
        self.hist_dict[('text content', 'readers')] = {
            'data': {rdr: len(tcids) for rdr, tcids in reader_tcids.items()},
            }
        self.hist_dict[('readings', 'readers')] = {
            'data': {rdr: len(rid_set) for rdr, rid_set in reader_rids.items()}
            }
        return

    def _make_histograms(self):
        """Plot all the histograms."""
        # Produce the histograms
        for (agged, agg_over), data in self.hist_dict.items():
            self._plot_hist(agged, agg_over, data['data'])
            label = '%s per %s' % (agged, agg_over)
            if not isinstance(data['data'], dict):
                arr = data['data']
                stat_dict = {'mean': arr.mean(), 'std': arr.std(),
                             'median': np.median(arr)}
                stat_str = ', '.join(['%s=%f' % (k, v)
                                      for k, v in stat_dict.items()])
                self.add_text(stat_str, style='Code', section='Plots')
                self.summary_dict[label.capitalize()] = stat_dict.copy()
        return

    def _make_text_summary(self):
        """Stash a text file summary of totals."""
        text_report_str = ''
        top_labels = ['Total readings', 'Content processed',
                      'Statements produced']
        for label in top_labels:
            text_str = '%s: %d\n' % (label, self.summary_dict[label])
            self.add_text(text_str, section='Summary Statistics')
            text_report_str += '%s: %d\n' % (label, self.summary_dict[label])

        for label, data in self.summary_dict.items():
            if label in top_labels:
                continue
            if isinstance(data, dict):
                text_report_str += '%s:\n' % label
                text_report_str += '\n'.join(['\t%s: %d' % (k, v)
                                              for k, v in data.items()])
                text_report_str += '\n'
            else:
                text_str = '%s: %d\n' % (label, data)
                self.add_text(text_str, section='Summary Statistics')
                text_report_str += text_str
        self.s3.put_object(Key=self.s3_prefix + 'summary.txt',
                           Body=text_report_str, Bucket=self.bucket_name)
        return

    def _make_job_line(self, key, value):
        """For job info section, produce one line."""
        self.add_text(key, section='Job Info', space=(1, 6))
        self.add_text(value, section='Job Info', style='Code')

    def report_statistics(self, reading_outputs, stmt_outputs, starts, ends):
        """Grab low-hanging-statistics and produce a pdf report.

        All the data generated is also stashed for future use in text and pickle
        files alongside the pdf summary on s3.

        Parameters
        ----------
        reading_outputs : list [ReadingData]
            A list of the ReadingData results output by `produce_readings`.
        stmt_outputs : list [StatementData]
            A list of the StatementData results output by `produce_statements`.
        starts : dict {<stage> : <datetime.datetime instance>}
            A dict of the start times of different parts of the job
            (e.g. statement production).
        ends : dict {<stage> : <datetime.datetime instance>}
            A dict of the end times of different parts of the job (like starts).
        """
        starts['stats'] = datetime.now()
        for k, end in ends.items():
            self._make_job_line(k + ' start', str(starts[k]))
            self._make_job_line(k + ' end', str(end))
            self._make_job_line(k + ' duration', str(end-starts[k]))

        self.summary_dict['Total readings'] = len(reading_outputs)
        reading_stmts = [(rd.reading_id, rd.tcid, rd.reader, rd.get_statements())
                         for rd in reading_outputs]
        self.summary_dict['Content processed'] = \
            len({t[1] for t in reading_stmts})
        self.summary_dict['Statements produced'] = len(stmt_outputs)
        self.s3.put_object(Key=self.s3_prefix + 'raw_tuples.pkl',
                           Body=pickle.dumps([t[:-1] + (len(t[-1]),)
                                              for t in reading_stmts]),
                           Bucket=self.bucket_name)

        readings_with_stmts = []
        readings_with_no_stmts = []
        for t in reading_stmts:
            if t[-1]:
                readings_with_stmts.append(t)
            else:
                readings_with_no_stmts.append(t)

        self.summary_dict['Readings with no statements'] = \
            len(readings_with_no_stmts)

        self._populate_hist_data(readings_with_stmts, readings_with_no_stmts)
        self._make_histograms()
        self._make_text_summary()
        ends['stats'] = datetime.now()

        fname = self.make_report()
        with open(fname, 'rb') as f:
            self.s3.put_object(Key=self.s3_prefix + fname, Body=f.read(),
                               Bucket=self.bucket_name)

        self._make_timing_report(starts, ends)
        self._stash_data()
        return
