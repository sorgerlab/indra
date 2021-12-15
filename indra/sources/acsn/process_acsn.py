import os
import gilda
import pandas as pd
from collections import Counter
from famplex import load_relations
from collections import defaultdict
from indra.statements.statements import stmts_from_json


def make_stmts(subj, obj, pmid, stmt_type):
    stmt = {'type': stmt_type,
            'subj': {'name': subj, 'db_refs': {'TEXT': 'kappa'}},
            'obj': {'name': obj, 'db_refs': {'TEXT': 'delta receptor'}},
            'obj_activity': 'activity',
            'evidence': [{'text_refs': {'PMID': pmid}}]}


HERE = os.path.join(os.path.realpath(os.path.dirname(__file__)))
acsn_df = pd.read_csv(os.path.join(HERE, 'ACSN2_binary_relations_between_proteins_with_PMID.txt'),
                      sep='\t')
acsn_gmt = os.path.join(HERE, 'ACSN2_HUGO_Correspondence.gmt')
famplex_df = pd.DataFrame(load_relations())
famplex_df.rename({0: 'ns_1', 1: 'subj', 2: 'relation', 3: 'ns_2', 4: 'obj'}, axis=1,
                  inplace=True)

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
            for gene in acsn_gmt_dict[ag]:
                if gene in list(famplex_df.subj):
                    idx = list(famplex_df.subj).index(gene)
                    grounded_gene = famplex_df.obj[idx]
                    grounded_db = 'FPLX'
                    agents[ag] = (grounded_gene, grounded_db)
                    break

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
        for t in text_ref.split(';'):
            stmt = {'type': stmt_type,
                    'subj': {'name': subj, 'db_refs': {subj_db_ref: subj}},
                    'obj': {'name': obj, 'db_refs': {obj_db_ref: obj}},
                    'obj_activity': 'activity',
                    'evidence': [{'text_refs': {'PMID': t}}]}
            stmts.append(stmt)

interaction_counts = Counter(acsn_df['INTERACTION_TYPE'])
