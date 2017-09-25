from util import *
import indra.tools.assemble_corpus as ac
from run_task5 import export_paths, report_paths

# This should be done for all cell lines
cell_lines = ['C32', 'RVH421']
for cell_line in cell_lines:
    scored_paths, models = pklload('task5_scored_paths')
    model = models[cell_line]
    stmts = pklload('pysb_stmts_%s' % cell_line)
    report_paths(scored_paths[cell_line], model, stmts, cell_line)
