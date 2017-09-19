import pickle
from copy import deepcopy
from indra.explanation.model_checker import ModelChecker
from process_data import read_rppa_data, get_task_1, get_antibody_agents
from assemble_models import prefixed_pkl
from assemble_pysb import contextualize_model


def get_agent_values(antibody_agents, values):
    agent_values = {}
    for ab, value in values.items():
        for agent in antibody_agents[ab]:
            agent_values[agent] = value
    return


def check_for_cell_line(cell_line, stmts, agents_to_observe):
    print('Cell line: %s\n=============' % cell_line)
    # Check for a single dose here since this part is not scored yet
    results = {}
    dose = 1.0
    for drug in stmts.keys():
        print('Drug: %s\n=============' % drug)
        stmts_condition = [stmt for stmt, _ in stmts[drug][dose]]
        mc = ModelChecker(model_cell_line, stmts_condition, agents_to_observe)
        path_results = mc.check_model(1, 5)
        results[drug] = path_results
    return path_results


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
    stmts_to_check = get_task_1(data)

    all_results = {}
    # Run model checking per cell line
    for cell_line in stmts_to_check.keys():
        model_cell_line = deepcopy(model)
        model_cell_line = contextualize_model(model_cell_line, cell_line)
        results = check_for_cell_line(cell_line, stmts_to_check[cell_line],
                                      agents_to_observe)
        all_results[cell_line] = results
    #agent_values = get_agent_values(antibody_agents, values)
