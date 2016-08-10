import os
import csv
from copy import deepcopy

class GroundingMapper(object):
    def __init__(self, gm):
        self.gm = gm

    def map_agents(self, stmts):
        # Make a copy of the stmts
        mapped_stmts = deepcopy(stmts)
        # Iterate over the statements
        for stmt in mapped_stmts:
            # Iterate over the agents
            for agent in stmt.agent_list():
                if agent is None or agent.db_refs.get('TEXT') is None:
                    continue
                agent_text = agent.db_refs.get('TEXT')
                # Look this string up in the grounding map
                agent_map_entry = self.gm.get(agent_text)
                # If not in the map, continue
                if agent_map_entry is None:
                    continue
                # Otherwise, update the agent's db_refs field
                agent.db_refs = agent_map_entry
        return mapped_stmts

def load_grounding_map(path):
    g_map = {}
    with open(path) as f:
        mapreader = csv.reader(f, delimiter='\t')
        for row in mapreader:
            key = row[0]
            db_refs = {'TEXT': key}
            for pair_ix in range(0, 2):
                col_ix = (pair_ix * 2) + 1
                db = row[col_ix]
                db_id = row[col_ix + 1]
                if db == '' or db == 'None' or db_id == '' or db_id == 'None':
                    continue
                else:
                    db_refs[db] = db_id
            if len(db_refs.keys()) > 1:
                g_map[key] = db_refs
    return g_map

default_grounding_map_path = os.path.join(os.path.dirname(__file__),
                                  '../resources/grounding_map.txt')
default_grounding_map = load_grounding_map(default_grounding_map_path)
gm = default_grounding_map

