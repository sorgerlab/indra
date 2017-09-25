from util import *
import indra.tools.assemble_corpus as ac
from run_task1 import export_paths, report_paths

cell_line = 'C32'
scored_paths, models = pklload('task1_scored_paths')
model = models[cell_line]
stmts = pklload('pysb_stmts_%s' % cell_line)

report_paths(scored_paths, model, stmts)
