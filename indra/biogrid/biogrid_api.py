import biogrid_client as bgc
from process import Publication

def process_query(agent_list, save_json_name = 'biogrid_output.json'):
    dict = bgc.get_json(agent_list)
    bgc.save_json(dict, save_json_name)
    return process_json(dict, agent_list)


def get_biogrid_interactors(dict):
    b_list = [dict['OFFICIAL_SYMBOL_A'], dict['OFFICIAL_SYMBOL_B']]
    return b_list
    

def get_subset(dict, agent_list):
    new_dict = {}
    for id in dict.keys():
        b_list = get_biogrid_interactors(dict[id])
        if set(b_list) == set(agent_list):
            new_dict[id] = js[id]
    return new_dict


def process_json(dict, agent_list):
    agent_subset = get_subset(dict, agent_list)
    publications = []
    for paper in agent_subset.keys():
        publications.append(Publication(agent_subset[paper]))
    return publications
        

            
