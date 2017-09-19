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

with open(prefixed_pkl('pysb_model'), 'rb') as fh:
    model = pickle.load(fh)

data = read_rppa_data()
antibody_agents = get_antibody_agents()
stmts_to_check = get_task_1(data)
for cell_line in stmts_to_check.keys():
    model_cell_line = deepcopy(model)
    model_cell_line = contextualize_model(model_cell_line, cell_line)
    for drug in stmts_to_check[cell_line].keys():
        for dose in stmts_to_check[cell_line][drug].keys():
            for stmt, values in stmts_to_check[cell_line][drug][dose]:
                agent_values = get_agent_values(antibody_agents, values)
                mc = ModelChecker(model_cell_line, [stmt], agent_values)
                path_results = mc.check_model(1, 5)
            break
        break
    break
