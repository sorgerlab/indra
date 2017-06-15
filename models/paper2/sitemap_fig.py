from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join as pjoin
from indra.tools import assemble_corpus as ac
from collections import OrderedDict, Counter
from indra.preassembler.sitemapper import SiteMapper, default_site_map
from indra.util import write_unicode_csv, read_unicode_csv

def get_incorrect_sites():
    outf = '../phase3_eval/output'

    prior_stmts = ac.load_statements(pjoin(outf, 'prior.pkl'))
    reach_stmts = ac.load_statements(pjoin(outf, 'phase3_stmts.pkl'))
    stmts = ac.map_grounding(prior_stmts + reach_stmts,
                             save=pjoin(outf, 'gmapped_stmts.pkl'))

    # Look for errors in database statements
    def get_inc_sites(stmts):
        sm = SiteMapper(default_site_map)
        valid, mapped = sm.map_sites(stmts)
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
        mapped = 1 if mm[1] else 0
        rows.append([gene, res, pos, freq, mapped])
    rows = sorted(rows, key=lambda x: x[3], reverse=True)
    write_unicode_csv('incorrect_sites.csv', rows)

def load_incorrect_sites():
    rows = read_unicode_csv('incorrect_sites.csv')
    return rows

if __name__ == '__main__':
    sites = load_incorrect_sites()

