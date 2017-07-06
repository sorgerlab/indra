from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import glob
import pickle
from os.path import join as pjoin
from indra.tools import assemble_corpus as ac
from collections import OrderedDict, Counter, namedtuple, defaultdict
from indra.preassembler.sitemapper import SiteMapper, default_site_map
from indra.util import write_unicode_csv, read_unicode_csv
from indra.util import plot_formatting as pf

import numpy as np
from matplotlib import pyplot as plt

pf.set_fig_params()


SiteInfo = namedtuple('SiteInfo', ['gene', 'res', 'pos', 'valid', 'mapped',
                                   'mapped_res', 'mapped_pos', 'explanation',
                                   'freq', 'source'])


def map_statements(stmts, source, outfile=None):
    """Tabulate valid, invalid, and mapped sites from a set of Statements."""
    # Look for errors in database statements
    sm = SiteMapper(default_site_map)
    valid_stmts, mapped_stmts = sm.map_sites(stmts)
    # Collect stats from SiteMapper itself
    sites = []
    for site_key, mapping in sm._cache.items():
        gene, res, pos = site_key
        freq = sm._sitecount[site_key]
        if mapping == 'VALID':
            valid, mapped, mapped_res, mapped_pos, explanation = \
                                                      (1, 0, None, None, None)
        else:
            valid = 0
            # Not mapped
            if mapping is None:
                mapped, mapped_res, mapped_pos, explanation = \
                                                    (0, None, None, None)
            # Mapped!
            else:
                mapped_res, mapped_pos, explanation = mapping
                mapped = 1 if mapped_pos else 0
        si = SiteInfo(gene, res, pos, valid, mapped, mapped_res, mapped_pos,
                      explanation, freq, source)
        sites.append(si)
    # Write to CSV file
    if outfile:
        header = [[field.upper() for field in si._asdict().keys()]]
        rows = header + replace_nones(sites)
        write_unicode_csv(outfile, rows)
    return sites


def map_agents(mod_agents_file, sm, source, save_csv=True):
    """Tabulate valid, invalid, and mapped sites from a set of Agents."""
    # Load the agents
    with open(mod_agents_file, 'rb') as f:
        mod_agents = pickle.load(f)
    print("Mapping %s" % mod_agents_file)
    sites = []
    for ag_ix, ag in enumerate(mod_agents):
        #if ag_ix % 1000 == 0:
        #    print('%d of %d' % (ag_ix, len(mod_agents)))
        invalid_sites = sm._check_agent_mod(ag, ag.mods, True, True, True)
        # Valid
        if not invalid_sites:
            valid, mapped, mapped_res, mapped_pos, explanation = \
                                                      (1, 0, None, None, None)
        else:
            assert len(invalid_sites) == 1
            mapping = invalid_sites[0][1]
            valid = 0
            # Not mapped
            if mapping is None:
                mapped, mapped_res, mapped_pos, explanation = \
                                                    (0, None, None, None)
            # Mapped!
            else:
                mapped_res, mapped_pos, explanation = mapping
                mapped = 1 if mapped_pos else 0
        si = SiteInfo(ag.name, ag.mods[0].residue, ag.mods[0].position, valid,
                      mapped, mapped_res, mapped_pos, explanation, None, source)
        sites.append(si)
    # Now that we've collected a list of all the sites, tabulate frequencies
    site_counter = Counter(sites)
    sites_with_freq = []
    for site, site_freq in site_counter.items():
        site_args = site._asdict()
        site_args.update({'freq': site_freq})
        si = SiteInfo(**site_args)
        sites_with_freq.append(si)
    # Write to CSV file
    if save_csv:
        header = [field.upper() for field in si._asdict().keys()]
        rows = header + replace_nones(sites)
        write_unicode_csv(agent_file.split('.')[0] + '.csv', rows)
    return sites_with_freq

# ------------------------------------------------------------

def load_incorrect_sites():
    rows = read_unicode_csv('incorrect_sites.csv')
    return [SiteInfo(*row) for row in rows]

def make_bar_plot(site_info, num_genes=120):
    # Build a dict based on gene name
    # Get counts summed across gene names
    gene_counts = defaultdict(lambda: 0)
    site_counts = defaultdict(list)
    for site in site_info:
        gene_counts[site.gene] += int(site.freq)
        site_counts[site.gene].append((int(site.freq), int(site.mapped),
                                      site.mapped_res, site.mapped_pos,
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
            for site_freq, mapped, mapped_res, mapped_pos, explanation \
                    in site_counts[gene]:
                if mapped and \
                        explanation.startswith('INFERRED_METHIONINE_CLEAVAGE'):
                    color = 'b'
                    handle_key = 'Methionine'
                elif mapped and not mapped_pos:
                    color = 'r'
                    handle_key = 'Curated as incorrect'
                elif mapped and \
                        explanation.startswith('INFERRED_MOUSE_SITE'):
                    color = 'c'
                    handle_key = 'Mouse'
                elif mapped and \
                        explanation.startswith('INFERRED_RAT_SITE'):
                    color = 'purple'
                    handle_key = 'Rat'
                elif mapped and \
                        explanation.startswith('INFERRED_ALTERNATIVE_ISOFORM'):
                    color = 'orange'
                    handle_key = 'Alternative isoform'
                elif mapped:
                    color = 'g'
                    handle_key = 'Manually mapped'
                else:
                    color = 'white'
                    handle_key = 'Unmapped'
                handle_dict[handle_key] = plt.bar(ix, site_freq, bottom=bottom,                                                   color=color,
                                                  linewidth=0.5, width=width)
                bottom += site_freq
        plt.xticks(ind + (width / 2.), [x[0] for x in gene_count_subset],
                   rotation='vertical')
        plt.ylabel('Stmts with invalid sites')
        plt.xlim((0, max(ind)+1))
        ax = plt.gca()
        pf.format_axis(ax)
        plt.subplots_adjust(**subplot_params)
        if do_legend:
            plt.legend(loc='upper right', handles=list(handle_dict.values()),
                       labels=list(handle_dict.keys()), fontsize=pf.fontsize,
                       frameon=False)
        plt.show()
    plot_sites(gene_counts[0:4], (0.23, 2),
               {'left': 0.24, 'right': 0.52, 'bottom': 0.31}, do_legend=False)
    plot_sites(gene_counts[4:num_genes], (11, 2),
               {'bottom': 0.31, 'left': 0.06, 'right':0.96})
    return gene_counts

# ---------------------


def replace_nones(rows):
    return [[cell if cell is not None else '' for cell in row]
             for row in rows]


def plot_pc_pe_mods(all_mods):
    dbs = Counter([row.source for row in all_mods])
    dbs = sorted([(k, v) for k, v in dbs.items()], key=lambda x: x[1],
                 reverse=True)
    plt.ion()
    width = 0.8
    ind = np.arange(len(dbs)) + (width / 2.)
    plt.figure(figsize=(2, 2), dpi=150)
    for db_ix, (db, db_freq) in enumerate(dbs):
        db_mods = [row for row in all_mods if row.source == db]
        valid = [row for row in db_mods if row.valid == 1]
        invalid = [row for row in db_mods if row.valid == 0]
        h_valid = plt.bar(db_ix, len(valid), width=width, color='g')
        h_invalid = plt.bar(db_ix, len(invalid), width=0.8, color='r',
                            bottom=len(valid))
    plt.xticks(ind, [db[0] for db in dbs])
    ax = plt.gca()
    pf.format_axis(ax)
    plt.show()

def plot_site_count_dist(sites, num_sites=240):
    # Plot site frequencies, colored by validity
    sites.sort(key=lambda s: s.freq, reverse=True)
    width = 0.8
    plt.ion()
    plt.figure(figsize=(11, 2), dpi=150)
    ind = np.arange(num_sites) + (width / 2.)
    for site_ix, site in enumerate(sites[:num_sites]):
        if site.valid:
            color = 'g'
        else:
            color = 'r'
        plt.bar(site_ix, site.freq, color=color)
    ax = plt.gca()
    pf.format_axis(ax)
    plt.show()

def print_stats(sites):
    # Tabulate some stats
    n = len(sites)
    n_val = len([s for s in sites if s.valid])
    n_inv = n - n_val
    n_map = len([s for s in sites if s.mapped])
    def pct(n, d):
        return 100 * n / float(d)
    f = np.sum([s.freq for s in sites])
    f_val = np.sum([s.freq for s in sites if s.valid])
    f_inv = f - f_val
    f_map = np.sum([s.freq for s in sites if s.mapped])
    print("Total sites: %d" % n)
    print("  Valid:   %d (%0.1f)" % (n_val, pct(n_val, n)))
    print("  Invalid: %d (%0.1f)" % (n_inv, pct(n_inv, n)))
    print("  Mapped:  %d (%0.1f)" % (n_map, pct(n_map, n)))
    print("%% Mapped:  %0.1f" % pct(n_map, n_inv))
    print()
    print("Total site occurrences: %d" % f)
    print("  Valid:   %d (%0.1f)" % (f_val, pct(f_val, f)))
    print("  Invalid: %d (%0.1f)" % (f_inv, pct(f_inv, f)))
    print("  Mapped:  %d (%0.1f)" % (f_map, pct(f_map, f)))
    print("Pct occurrences mapped: %0.1f" % pct(f_map, f_inv))
    print()
    # Sample 100 invalid-unmapped (by unique sites)
    # Sample 100 invalid-mapped (by unique sites)

if __name__ == '__main__':
    outf = '../phase3_eval/output'
    prior_stmts = ac.load_statements(pjoin(outf, 'prior.pkl'))
    site_info = map_statements(prior_stmts, source='prior',
                               outfile='prior_sites.csv')

    #reach_stmts = ac.load_statements(pjoin(outf, 'phase3_stmts.pkl'))
    #stmts = prior_stmts
    #stmts = reach_stmts
    #stmts = ac.map_grounding(stmts, save=pjoin(outf, 'gmapped_stmts.pkl'))
    #stmts = ac.load_statements(pjoin(outf, 'gmapped_stmts.pkl'))

    sys.exit()
    """
    valid, sites, sm = get_incorrect_sites(do_methionine_offset=True,
                                 do_orthology_mapping=True,
                                 do_isoform_mapping=True)
    with open('sm.pkl', 'wb') as f:
        pickle.dump((sm._cache, sm._sitecount), f)
    plot_site_count_dist(sm)
    gene_counts = make_bar_plot(sites)

    """
    # This script does two things:
    # 1) Plots stats on invalid sites from databases
    #    - showing their frequency
    #       - per site
    #       - per reaction
    # 2) Showing the fraction of the invalid sites in DBs that are mapped
    #    - per site
    #    - per reaction
    # 3) Showing accuracy:
    #    - that the mapped sites are likely legit
    #    - and that the unmapped sites are likely errors

    sm = SiteMapper(default_site_map)
    with open('smcache.pkl', 'rb') as f:
        (sm._cache, sm._sitecount) = pickle.load(f)

    # Load the agent files
    agent_files = ['pc_pid_modified_agents.pkl',
                   'pc_psp_modified_agents.pkl',
                   'pc_reactome_modified_agents.pkl']
    # For each set of mods
    all_sites = []
    for agent_file in agent_files:
        db_name = agent_file.split('_')[1]
        sites = map_agents(agent_file, sm, db_name)
        all_sites += sites
        print("Stats for %s -------------" % db_name)
        print_stats(sites)
    header = [field.upper() for field in all_sites[0]._asdict().keys()]
    rows = header + replace_nones(all_sites)
    write_unicode_csv('pc_all_pe_mods.csv', rows)
    plot_pc_pe_mods(all_sites)
    #with open('smcache.pkl', 'wb') as f:
    #    pickle.dump((sm._cache, sm._sitecount), f)



