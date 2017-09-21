import pickle
from copy import deepcopy
import itertools
from indra.explanation.model_checker import ModelChecker
from process_data import read_rppa_data, get_task_1, get_antibody_agents
from assemble_models import prefixed_pkl
from assemble_pysb import contextualize_model


def get_agent_values(antibody_agents, values):
    agent_values = {}
    for ab, value in values.items():
        for agent in antibody_agents[ab]:
            agent_values[agent] = value
    return agent_values


flatten = lambda x: list(itertools.chain.from_iterable(x))

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
            mc = ModelChecker(model_cell_line, stmts_condition,
                              agents_to_observe)
            path_results = mc.check_model(10, 5)
            # Get a dict of measured values by INDRA Agents for this condition
            agent_values = get_agent_values(antibody_agents, values_condition)
            # Get a single list of paths for this condition
            # The length of this will be the max path number times the number
            # of Statements being checked
            paths = flatten([pr[1].paths for pr in path_results])
            # Run score paths with
            # - the list of paths the model checker produced
            # - the values of observed agents in the current condition
            # - loss of function mode since these are inhibitory drugs
            # - sigma corresponding to the log2 normalized measurements
            scored_paths_condition = mc.score_paths(paths, agent_values,
                                                    loss_of_function=True,
                                                    sigma=0.5)
            scored_paths[cell_line][drug] = scored_paths_condition

