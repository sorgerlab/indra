import os
import csv


class GroundingMapper(object):
    def __init__(self, gm):
        self.gm = gm

    def map_agents(self, stmts):
        pass

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
                                  '../resources/grounding_map.tsv')
default_grounding_map = load_grounding_map(default_grounding_map_path)
gm = default_grounding_map

