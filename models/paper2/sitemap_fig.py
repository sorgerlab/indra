from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join as pjoin
from indra.tools import assemble_corpus as ac
from collections import OrderedDict, Counter, namedtuple, defaultdict
from indra.preassembler.sitemapper import SiteMapper, default_site_map
from indra.util import write_unicode_csv, read_unicode_csv
from indra.util import plot_formatting as pf

import numpy as np
from matplotlib import pyplot as plt

pf.set_fig_params()


SiteInfo = namedtuple('SiteInfo', ['gene', 'res', 'pos', 'freq', 'mapped',
                                   'mapped_res', 'mapped_pos', 'explanation'])

def get_incorrect_sites(do_methionine_offset=False, do_orthology_mapping=False):
    outf = '../phase3_eval/output'

    prior_stmts = ac.load_statements(pjoin(outf, 'prior.pkl'))
    #reach_stmts = ac.load_statements(pjoin(outf, 'phase3_stmts.pkl'))
    stmts = prior_stmts
    stmts = ac.map_grounding(stmts, save=pjoin(outf, 'gmapped_stmts.pkl'))

    # Look for errors in database statements
    def get_inc_sites(stmts):
        sm = SiteMapper(default_site_map)
        valid, mapped = sm.map_sites(stmts,
                                     do_methionine_offset=do_methionine_offset,
                                     do_orthology_mapping=do_orthology_mapping)
        # Collect stats on most frequently occurring site errors
        mapped_sites = []
        unmapped_sites = []
        for ms in mapped:
            for mm in ms.mapped_mods:
                if mm[1]:
                    mapped_sites.append(mm)
                else:
                    unmapped_sites.append(mm)
        unmapped_count = Counter(unmapped_sites)
        unmapped_sites = sorted([(k, v) for k, v in unmapped_count.items()],
                                 key=lambda x: x[1], reverse=True)
        mapped_count = Counter(mapped_sites)
        mapped_sites = sorted([(k, v) for k, v in mapped_count.items()],
                              key=lambda x: x[1], reverse=True)
        return (unmapped_sites, mapped_sites)

    (db_inc, db_map) = get_inc_sites(stmts)

    rows = []
    for mm, freq in db_inc + db_map:
        gene, res, pos = mm[0]
        if mm[1]:
            mapped = 1
            (mapped_res, mapped_pos, explanation) = mm[1]
        else:
            mapped = 0
            mapped_res = mapped_pos = explanation = ''
        rows.append([gene, res, pos, freq, mapped, mapped_res, mapped_pos,
                     explanation])
    rows = sorted(rows, key=lambda x: x[3], reverse=True)
    write_unicode_csv('incorrect_sites.csv', rows)
    return [SiteInfo(*row) for row in rows]

def load_incorrect_sites():
    rows = read_unicode_csv('incorrect_sites.csv')
    return [SiteInfo(*row) for row in rows]

def make_bar_plot(site_info, num_genes=60):
    # Build a dict based on gene name
    # Get counts summed across gene names
    gene_counts = defaultdict(lambda: 0)
    site_counts = defaultdict(list)
    for site in site_info:
        gene_counts[site.gene] += int(site.freq)
        site_counts[site.gene].append((int(site.freq), int(site.mapped),
                                       site.explanation))
    # Sort the individual site counts by frequency
    for gene, freq_list in site_counts.items():
        site_counts[gene] = sorted(freq_list, key=lambda x: x[0], reverse=True)
    gene_counts = sorted([(k, v) for k, v in gene_counts.items()],
                          key=lambda x: x[1], reverse=True)

    plt.ion()
    def plot_sites(gene_count_subset, figsize, subplot_params, do_legend=True):
        ind = np.array(range(len(gene_count_subset)))
        plt.figure(figsize=figsize, dpi=150)
        width = 0.8
        handle_dict = {}
        for ix, (gene, freq) in enumerate(gene_count_subset):
            # Plot the stacked bars
            bottom = 0
            for site_freq, mapped, explanation in site_counts[gene]:
                if mapped and \
                        explanation.startswith('INFERRED_METHIONINE_CLEAVAGE'):
                    color = 'b'
                    handle_key = 'Methione'
                elif mapped and \
                        explanation.startswith('INFERRED_MOUSE_SITE'):
                    color = 'r'
                    handle_key = 'Mouse'
                elif mapped and \
                        explanation.startswith('INFERRED_RAT_SITE'):
                    color = 'purple'
                    handle_key = 'Rat'
                elif mapped:
                    color = 'g'
                    handle_key = 'Manual'
                else:
                    color = 'white'
                    handle_key = 'Unmapped'
                handle_dict[handle_key] = plt.bar(ix, site_freq, bottom=bottom,                                                   color=color,
                                                  linewidth=0.5, width=width)
                bottom += site_freq
        plt.xticks(ind + (width / 2.), [x[0] for x in gene_count_subset],
                   rotation='vertical')
        plt.ylabel('Num. invalid sites')
        plt.xlim((0, max(ind)+1))
        ax = plt.gca()
        pf.format_axis(ax)
        plt.subplots_adjust(**subplot_params)
        if do_legend:
            plt.legend(loc='upper right', handles=list(handle_dict.values()),
                       labels=list(handle_dict.keys()), fontsize=pf.fontsize,
                       frameon=False)
        plt.show()
    plot_sites(gene_counts[0:4], (1, 2), {'bottom': 0.31}, do_legend=False)
    plot_sites(gene_counts[4:num_genes], (7, 2),
               {'bottom': 0.31, 'left': 0.06, 'right':0.96})
    return gene_counts

if __name__ == '__main__':
    #sites = load_incorrect_sites()
    sites = get_incorrect_sites(do_methionine_offset=True,
                                do_orthology_mapping=True)
    gene_counts = make_bar_plot(sites)
