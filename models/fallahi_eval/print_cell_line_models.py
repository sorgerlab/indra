from util import *
from process_data import *
from assemble_pysb import print_initial_conditions

data = read_rppa_data()
gene_names = get_gene_names(data)

fname = os.path.join(based, basen + '_cell_line_parameters.tsv')

print_initial_conditions(cell_lines, gene_names, fname)
