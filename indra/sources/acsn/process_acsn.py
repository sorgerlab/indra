import os
import pandas as pd
from collections import defaultdict, Counter

from indra.ontology.bio import bio_ontology


def get_famplex_lookup():
    fplx_lookup = {}
    bio_ontology.initialize()
    for node in bio_ontology.nodes:
        ns, id = bio_ontology.get_ns_id(node)
        if ns == 'FPLX':
            children = bio_ontology.get_children(ns, id)
            hgnc_children = [bio_ontology.get_name(*c)
                             for c in children if c[0] == 'HGNC']
            fplx_lookup[tuple(sorted(hgnc_children))] = id

fplx_lookup = get_famplex_lookup()

HERE = os.path.join(os.path.realpath(os.path.dirname(__file__)))
acsn_df = pd.read_csv(os.path.join(HERE, 'ACSN2_binary_relations_between_proteins_with_PMID.txt'),
                      sep='\t')
acsn_gmt = os.path.join(HERE, 'ACSN2_HUGO_Correspondence.gmt')

# Convert the GMT file into a dictionary
acsn_gmt_dict = defaultdict(set)
with open(acsn_gmt, 'r') as fh:
    for line in fh:
        line = line.rstrip('\n')
        line = line.rstrip('\t')
        interactors = line.split('\t')
        for ag2 in interactors[2:]:
            acsn_gmt_dict[(interactors[0])].add(ag2)

stmts = []
count = 0
unmapped_genes = []
# Create INDRA statements
for v in acsn_df.values:
    int_1, int_2 = v[0], v[2]
    stmt_type = v[1]
    text_ref = v[3]
    agents = {int_1: 'NA',
              int_2: 'NA'}

    # Grounding the agents
    for ag in [int_1, int_2]:
        if ag in acsn_gmt_dict and len(acsn_gmt_dict[ag]) == 1:
            grounded_gene = next(iter(acsn_gmt_dict[ag]))
            grounded_db = 'HGNC'
            agents[ag] = (grounded_gene, grounded_db)

        elif ag in acsn_gmt_dict and len(acsn_gmt_dict[ag]) > 1:
            fplx_id = fplx_lookup.get(tuple(sorted(acsn_gmt_dict[ag])))
            if fplx_id:
                agents[ag] = (fplx_id, 'FPLX')

    # Check for the unmapped agents and exclude that statement
    if agents[int_1] == 'NA' or agents[int_2] == 'NA':
        count += 1
        if agents[int_1] == 'NA':
            unmapped_genes.append(int_1)
        if agents[int_2] == 'NA':
            unmapped_genes.append(int_2)
        continue

    # create and append the statement into statements list
    subj = agents[int_1][0]
    subj_db_ref = agents[int_1][1]
    obj = agents[int_2][0]
    obj_db_ref = agents[int_2][1]

    # Check for Nan values in text_ref
    if str(text_ref) == 'nan':
        stmt = {'type': stmt_type,
                'subj': {'name': subj, 'db_refs': {subj_db_ref: subj}},
                'obj': {'name': obj, 'db_refs': {obj_db_ref: obj}},
                'obj_activity': 'activity',
                'evidence': [{'text_refs': {}}]}
        stmts.append(stmt)
    else:
        for id in text_ref.split(';'):
            stmt = {'type': stmt_type,
                    'subj': {'name': subj, 'db_refs': {subj_db_ref: subj}},
                    'obj': {'name': obj, 'db_refs': {obj_db_ref: obj}},
                    'obj_activity': 'activity',
                    'evidence': [{'text_refs': {'PMID': id}}]}
            stmts.append(stmt)

interaction_counts = Counter(acsn_df['INTERACTION_TYPE'])
