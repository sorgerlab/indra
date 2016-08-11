import os
import csv
from copy import deepcopy
from indra.databases import uniprot_client

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

    def rename_agents(self, stmts):
        # Make a copy of the stmts
        mapped_stmts = deepcopy(stmts)
        # Iterate over the statements
        for stmt_ix, stmt in enumerate(mapped_stmts):
            # Iterate over the agents
            for agent in stmt.agent_list():
                if agent is None:
                    continue
                old_name = agent.name
                # If there's an INDRA ID, prefer that for the name
                if agent.db_refs.get('INDRA'):
                    agent.name = agent.db_refs.get('INDRA')
                # Take a HGNC name from Uniprot next
                elif agent.db_refs.get('UP'):
                    # Try for the HGNC name
                    hgnc_name = uniprot_client.get_hgnc_name(
                                                    agent.db_refs.get('UP'))
                    if hgnc_name is not None:
                        agent.name = hgnc_name
                        continue
                    # Fall back on the Uniprot gene name
                    up_gene_name = uniprot_client.get_gene_name(
                                                    agent.db_refs.get('UP'))
                    if up_gene_name is not None:
                        agent.name = up_gene_name
                        continue
                    # Take the text string
                    if agent.db_refs.get('TEXT'):
                        agent.name = agent.db_refs.get('TEXT')
                    # If this fails, then we continue with no change
                # Fall back to the text string
                elif agent.db_refs.get('TEXT'):
                    agent.name = agent.db_refs.get('TEXT')
                if old_name != agent.name:
                    print "Map %d of %d: %s --> %s" % \
                                (stmt_ix+1, len(stmts), old_name, agent.name)
        return mapped_stmts

# TODO: handle the cases when there is more than one entry for the same
# key (e.g., ROS, ER)
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

