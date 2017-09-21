import pickle
from copy import deepcopy
import itertools
import indra.tools.assemble_corpus as ac
from indra.explanation.model_checker import ModelChecker, _stmt_from_rule
from indra.explanation.reporting import *
from process_data import read_rppa_data, get_task_1, get_antibody_agents
from assemble_pysb import contextualize_model, prefixed_pkl


def get_agent_values(antibody_agents, values):
    agent_values = {}
    for ab, value in values.items():
        for agent in antibody_agents[ab]:
            agent_values[agent] = value
    return agent_values


def get_global_mc(model, stmts_to_check, agents_to_observe):
    all_stmts_condition = []
    for cell_line in stmts_to_check.keys():
        for drug in stmts_to_check[cell_line].keys():
            stmts_condition, _ = stmts_to_check[cell_line][drug][1.0]
            all_stmts_condition += stmts_condition
    mc = ModelChecker(model, all_stmts_condition, agents_to_observe)
    mc.prune_influence_map()
    return mc

def export_paths(scored_paths, model):
    stmts = ac.load_statements(prefixed_pkl('pysb_stmts'))
    paths = get_paths(scored_paths, model, stmts)

flatten = lambda x: list(itertools.chain.from_iterable(x))

"""
def group_scored_paths(scored_paths, mc):
    groups = set()
    for path, score in scored_paths:
        for ix, rule in enumerate(path):
            # The last rule should be an observable, which won't have a
            # statement associated with it (for now)
            if ix == len(path) - 1:
                rule_stmt = path[-2]
                obj = rule_stmt.agent_list()[1]
                gene_name = obj.db_refs.get('HGNC')
            else:
                rule_stmt = _stmt_from_rule(mc.model, rule, mc.statements)
                if not rule_stmt:
                    continue
                subj = rule_stmt.agent_list()[0]
                gene_name = subj.db_refs.get('HGNC')
            if gene_name:
                grouped_path.append(gene_name)
        # At end, get the object of the last rule
"""

if __name__ == '__main__':
    # Some basic setup
    with open(prefixed_pkl('pysb_model'), 'rb') as fh:
        model = pickle.load(fh)
    data = read_rppa_data()
    antibody_agents = get_antibody_agents()
    agents_to_observe = []
    for agents in antibody_agents.values():
        agents_to_observe += agents

    # Get all the data for Task 1
    stmts_to_check = get_task_1(data, False)
    global_mc = None
    dose = 1.0
    scored_paths = {}
    # Run model checking per cell line
    for cell_line in stmts_to_check.keys():
        print('Cell line: %s\n=============' % cell_line)
        model_cell_line = deepcopy(model)
        model_cell_line = contextualize_model(model_cell_line, cell_line)
        scored_paths[cell_line] = {}
        for drug in stmts_to_check[cell_line].keys():
            print('Drug: %s\n=============' % drug)
            # Get all the statements and values for the condition
            stmts_condition, values_condition = \
                stmts_to_check[cell_line][drug][dose]
            # Make a Model Checker with the
            # - model contextualized to the cell line
            # - the statements for the given condition
            # - agents for which observables need to be made
            #mc = ModelChecker(model_cell_line, stmts_condition,
            #                  agents_to_observe)
            if not global_mc:
                global_mc = get_global_mc(model_cell_line, stmts_to_check,
                                          agents_to_observe)
            path_results = []
            for stmt in stmts_condition:
                pr = global_mc.check_statement(stmt, 20, 10)
                path_results.append(pr)
            #path_results = mc.check_model(2, 5)
            # Get a dict of measured values by INDRA Agents for this condition
            agent_values = get_agent_values(antibody_agents, values_condition)
            # Get a single list of paths for this condition
            # The length of this will be the max path number times the number
            # of Statements being checked
            paths = flatten([pr.paths for pr in path_results])
            # Run score paths with
            # - the list of paths the model checker produced
            # - the values of observed agents in the current condition
            # - loss of function mode since these are inhibitory drugs
            # - sigma corresponding to the log2 normalized measurements
            scored_paths_condition = global_mc.score_paths(paths, agent_values,
                                                    loss_of_function=True,
                                                    sigma=0.5)
            scored_paths[cell_line][drug] = scored_paths_condition



