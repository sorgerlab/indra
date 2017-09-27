from util import *
from run_task5 import report_paths
import group_paths as gp

stmts = pklload('pysb_stmts')

for pkl_name in ('task5_scored_paths', 'task5_scored_paths_inverse'):
    if pkl_name.endswith('_inverse'):
        print("\n\n -- INVERSE -------------------------------------\n\n")
    scored_path_dict, models = pklload(pkl_name)
    for cell_line in ['C32', 'RVH421']:
        model = models[cell_line]
        report_paths(scored_path_dict[cell_line], model, stmts, cell_line)

        for drug_name, scored_paths in scored_path_dict[cell_line].items():
            print("Grouping paths for %s, %s" % (cell_line, drug_name))
            gp.print_top_group_scores(scored_paths, model, stmts)


