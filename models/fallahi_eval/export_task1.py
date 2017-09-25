from util import *
import indra.tools.assemble_corpus as ac
from run_task1 import export_paths, report_paths
import group_paths as gp

"""
for cell_line in ['C32', 'LOXIMVI', 'MMACSF', 'MZ7MEL', 'RVH421']:
    model = models[cell_line]
    stmts = pklload('pysb_stmts_%s' % cell_line)
    report_paths(scored_path_dict[cell_line], model, stmts, cell_line)
"""


def path_group_scores(scored_path_dict, models):
    for cell_line in ['C32']: #, 'LOXIMVI', 'MMACSF', 'MZ7MEL', 'RVH421']:
        model = models[cell_line]
        stmts = pklload('pysb_stmts')
        for drug_name, scored_paths in scored_path_dict[cell_line].items():
            print("Grouping paths for %s, %s" % (cell_line, drug_name))
            groups, path_details = gp.group_scored_paths(scored_paths,
                                                         model, stmts)
            all_scores = gp.all_scores(path_details)
            np = gp.num_paths_by_group(path_details)
            ts = gp.top_scores(path_details)
            if cell_line == 'C32' and drug_name == 'Vemurafenib':
                gp.print_list(ts)
                for path, score in ts:
                    path_str = ' -> '.join([p[0] for p in path])
                    print('%s : score %s' % (path_str, score))

if __name__ == '__main__':
    scored_path_dict, models = pklload('task1_scored_paths')
    path_group_scores(scored_path_dict, models)

    scored_path_dict, models = pklload('task1_scored_paths_inverse')
    path_group_scores(scored_path_dict, models)

