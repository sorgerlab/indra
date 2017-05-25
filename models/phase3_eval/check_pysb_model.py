import pickle
import itertools
from indra.util import write_unicode_csv
from indra.assemblers import PysbAssembler
from indra.tools.model_checker import ModelChecker
import indra.tools.assemble_corpus as ac
import process_data
import make_stmts_for_checking as make_stmts
from assemble_pysb import set_context, add_observables

print("Processing data")

data = process_data.read_data(process_data.data_file)
data_genes = process_data.get_all_gene_names(data)
ab_map = process_data.get_antibody_map(data)

print('Loading data statements.')
data_stmts, data_values = make_stmts.run(dec_thresh=0.8, inc_thresh=1.2)
all_data_stmts = [values.values() for values in data_stmts.values()]
all_data_stmts = itertools.chain.from_iterable(all_data_stmts)
all_data_stmts = list(itertools.chain.from_iterable(all_data_stmts))


agent_obs = list(itertools.chain.from_iterable(ab_map.values()))
# Here we need to cross-reference the antbody map with the data values
agent_data = {}
for drug_name, values in data_values.items():
    agent_data[drug_name] = {}
    for ab_name, value in values.items():
        agents = ab_map[ab_name]
        for agent in agents:
            agent_data[drug_name][agent] = value

"""
base_stmts = ac.load_statements('output/korkut_model_pysb_no_evidence.pkl')

# Merge the sources of statements
# stmts = manual_stmts + base_stmts
stmts = base_stmts
#stmts = manual_stmts

# Assemble model
pa = PysbAssembler()
pa.add_statements(stmts)
model = pa.make_model()

with open('korkut_pysb.pkl', 'wb') as f:
    print("Pickling PySB model")
    pickle.dump(pa.model, f)
"""

with open('korkut_pysb.pkl', 'rb') as f:
    print("Unpickling PySB model")
    model = pickle.load(f)

# Preprocess and assemble the pysb model
#model = assemble_pysb(combined_stmts, data_genes, '')

mc = ModelChecker(model, all_data_stmts, agent_obs)

# Iterate over each drug/ab statement subset
results = []
for drug_name, ab_dict in data_stmts.items():
    agent_values = agent_data[drug_name]
    for ab, stmt_list in ab_dict.items():
        value = data_values[drug_name][ab]
        # For each subset, check statements; if any of them checks out, we're
        # good and can move on to the next group
        print("-- Checking the effect of %s on %s --" % (drug_name, ab))
        relation = 'positive' if value > 1 else 'negative'
        path_found = 0
        paths = []
        for stmt in stmt_list:
            print("Checking: %s" % stmt)
            result = mc.check_statement(stmt, max_paths=5, max_path_length=5)
            if result:
                paths += result[1]
                print("Path found, skipping rest")
                path_found = 1
                # Score paths here TODO
                #scored_result = mc.score_paths(result[0], agent_values)
                break
            else:
                print("No path found")

        results.append((drug_name, ab, relation, value, path_found, paths))
#write_unicode_csv('model_check_results.csv', results)
