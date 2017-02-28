import pickle
from indra.tools.model_checker import ModelChecker
from manual_stmts import stmts as manual_stmts
from assemble_pysb import assemble_pysb
import process_data
from indra.util import write_unicode_csv

data = process_data.read_data(process_data.data_file)
data_genes = process_data.get_all_gene_names(data)

with open('data_stmts.pkl', 'rb') as f:
    print('Loading data statements.')
    data_stmts, data_values = pickle.load(f)

#with open('output/korkut_model_pysb.pkl', 'rb') as f:
with open('korkut_stmts_no_ev.pkl', 'rb') as f:
    print('Loading korkut_model_pysb tatements.')
    base_stmts = pickle.load(f)

# Merge the sources of statements
#combined_stmts = manual_stmts + base_stmts
combined_stmts = base_stmts

# Preprocess and assemble the pysb model
model = assemble_pysb(combined_stmts, data_genes, '')

mc = ModelChecker(model)

# Iterate over each drug/ab statement subset
results = []
for drug_name, ab_dict in data_stmts.items():
    if drug_name not in ('901',):
        continue
    for ab, stmt_list in ab_dict.items():
        # For each subset, check statements; if any of them checks out, we're
        # good and can move on to the next group
        print("-- Checking the effect of %s on %s --" % (drug_name, ab))
        value = data_values[drug_name][ab]
        relation = 'positive' if value > 1 else 'negative'
        path_found = 0
        for stmt in stmt_list:
            print("Checking: %s" % stmt)
            result = mc.check_statement(stmt)
            if result:
                print("Path found, skipping rest")
                path_found = 1
                break
            else:
                print("No path found")

        results.append((drug_name, ab, relation, value, path_found))
write_unicode_csv('model_check_results.csv')
