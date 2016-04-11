import csv
from indra.bel import bel_api
from indra.biopax import biopax_api as ba
from indra.preassembler import Preassembler, render_stmt_graph
from indra.preassembler.hierarchy_manager import entity_hierarchy as eh, \
                                                 modification_hierarchy as mh
import pickle
from indra.preassembler.sitemapper import default_mapper as sm
from indra.statements import *
from pysb import kappa
from indra.pysb_assembler import PysbAssembler

gene_list = []

# STEP 0: Get gene list
# Get gene list from ras_pathway_proteins.csv
with open('data/ras_pathway_proteins.csv') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        gene_list.append(row[0].strip())

# STEP 1: Get statements from BEL
bel_statements = []
for gene in gene_list:
    print "Getting statements for gene", gene
    bel_proc = bel_api.process_ndex_neighborhood([gene])
    bel_statements += bel_proc.statements

with open('ras_genes_bel_stmts.pkl', 'w') as f:
    pickle.dump(bel_statements, f)

# STEP 2: Get Biopax statements
# Check for cached file before querying Pathway Commons Web API
biopax_ras_owl_path = 'ras_genes_pc_pathsbetween.owl'
if os.path.isfile(biopax_ras_owl_path):
    print "Loading Biopax from OWL file", biopax_ras_owl_path
    bp = ba.process_owl(biopax_ras_owl_path)
# save to owl file
else:
    print "Querying Pathway Commons for pathsbetween genes"
    bp = ba.process_pc_pathsbetween(gene_list)
    bp.save_model(biopax_ras_owl_path)

bp.get_phosphorylation()
bp.get_dephosphorylation()
bp.get_acetylation()
bp.get_palmitoylation()
bp.get_glycosylation()
bp.get_activity_modification()

stmts = bp.statements + bel_statements
with open('ras_genes_stmts.pkl', 'w') as f:
    pickle.dump(pa.unique_stmts, f)
"""

# STEP 3
with open('ras_genes_stmts.pkl') as f:
    stmts = pickle.load(f)

stmts = [s for s in stmts if isinstance(s, Modification)]

pa1 = Preassembler(eh, mh, stmts)
print "Combining duplicates"
pa1.combine_duplicates()

print "Mapping sites"
(valid, mapped) = sm.map_sites(pa1.unique_stmts)

mapped_stmts = valid + [m.mapped_stmt for m in mapped]

pa2 = Preassembler(eh, mh, mapped_stmts)
print "Combining duplicates again"
pa2.combine_duplicates()

import sys; sys.exit()

print "Combining related"
pa2.combine_related()

results = {}
results['raw'] = stmts
results['duplicates1'] = pa1.unique_stmts
results['valid'] = valid
results['mapped'] = mapped
results['mapped_stmts'] = mapped_stmts
results['duplicates2'] = pa2.unique_stmts
results['related2'] = pa2.related_stmts

with open('ras_genes_results.pkl', 'w') as f:
    pickle.dump(results, f)


# STEP 4
with open('ras_genes_results.pkl') as f:
    results = pickle.load(f)

sublist = ['RAF1', 'MAP2K1', 'MAPK1', 'KSR1', 'DUSP1', 'KRAS', 'AKT1', 'PDPK1']

subgraph = []
for s in results['related2']:
    in_sub_list = [agent is None or agent.name in sublist
                   for agent in s.agent_list()]
    if all(in_sub_list):
        subgraph.append(s)

with open('ras_genes_subgraph.pkl', 'w') as f:
    pickle.dump(subgraph, f)

# STEP 5 Kappa contact map
with open('ras_genes_subgraph.pkl') as f:
    stmts = pickle.load(f)
stmts = [s for s in stmts \
         if not (s.enz.name == 'MAPK1' and s.sub.name == 'RAF1')]
stmts = [s for s in stmts \
         if not (s.enz.name == 'PDPK1' and s.sub.name == 'AKT1')]


pya = PysbAssembler(policies='interactions_only')
pya.add_statements(stmts)
pya.make_model()
pya.model.name = 'ras_subgraph'
g = kappa.contact_map(pya.model)
g.draw('ras_genes_subgraph_cm.pdf', prog='dot')

# STEP 6
# Load the results from reading 100 papers
with open('pmc_papers/preassembler.pkl') as f:
    braf_stmts = pickle.load(f)

pa3 = Preassembler(eh, mh, braf_stmts)
pa3.combine_duplicates()

# Filter statements based on whether they were grounded to a gene/protein
braf_grounded = []
for s in pa3.unique_stmts:
    hgnc_grounded = [agent is None or
                     agent.db_refs.get('HGNC', None) is not None or
                     agent.db_refs.get('UP', None) is not None
                     for agent in s.agent_list()]
    if all(hgnc_grounded):
        braf_grounded.append(s)

# Filter statements based on whether the entities are in the McCormick list
braf_filtered = []
for s in braf_grounded:
    in_ras_list = [agent is None or agent.name in gene_list
                   for agent in s.agent_list()]
    if all(in_ras_list):
        braf_filtered.append(s)

#from indra.sbgn_assembler import SBGNAssembler
#sbgn = SBGNAssembler(rel)
#sbgn.make_sbgn()
"""
