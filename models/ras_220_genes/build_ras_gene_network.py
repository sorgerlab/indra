from indra.tools.gene_network import GeneNetwork, grounding_filter
import csv

# STEP 0: Get gene list
gene_list = []
# Get gene list from ras_pathway_proteins.csv
with open('../../data/ras_pathway_proteins.csv') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        gene_list.append(row[0].strip())

gn = GeneNetwork(gene_list, 'ras_genes')
stmts = gn.get_statements(filter=True)
grounded_stmts = grounding_filter(stmts)
results = gn.run_preassembly(grounded_stmts)

