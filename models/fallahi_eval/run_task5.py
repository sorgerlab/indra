import pickle
from copy import deepcopy
import itertools
from indra.explanation.model_checker import ModelChecker, _stmt_from_rule
from util import pklload, pkldump
from process_data import read_rppa_data, get_task_5, get_antibody_agents
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

flatten = lambda x: list(itertools.chain.from_iterable(x))

if __name__ == '__main__':
    INVERSE = False
    if INVERSE:
        print('Running Task 5 in INVERSE mode as control')
        print('=========================================')
    data = read_rppa_data()
    antibody_agents = get_antibody_agents()
    agents_to_observe = []
    for agents in antibody_agents.values():
        agents_to_observe += agents

    # Get all the data for Task 5
    stmts_to_check = get_task_5(data, INVERSE)
    dose = 1.0
    scored_paths = {}
    models = {}
    global_mc = None
    # Run model checking per cell line
    for cell_line in stmts_to_check.keys():
        print('Cell line: %s\n=============' % cell_line)
        model = pklload('pysb_model_%s' % cell_line)
        scored_paths[cell_line] = {}
        # Make a Model Checker with the
        # - model contextualized to the cell line
        # - the statements for the given condition
        # - agents for which observables need to be made
        global_mc = get_global_mc(model, stmts_to_check, agents_to_observe)
        for drug in stmts_to_check[cell_line].keys():
            print('Drug: %s\n=============' % drug)
            # Get all the statements and values for the condition
            stmts_condition, values_condition = \
                stmts_to_check[cell_line][drug][dose]
            path_results = []
            for stmt in stmts_condition:
                pr = global_mc.check_statement(stmt, 1, 10)
                path_results.append(pr)
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
        models[cell_line] = global_mc.model

# Dump results in standard folder
fname = 'task5_scored_paths'
if INVERSE:
    fname += '_inverse'
pkldump(fname, (scored_paths, models))
