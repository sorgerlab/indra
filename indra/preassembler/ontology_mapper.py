import os

class OntologyMapper(object):
    """A class to map between ontologies in grounded arguments of Statements.

    Parameters
    ----------
    statements : list[indra.statement.Statement]
        A list of INDRA Statements to map
    mappings : Optional[list[tuple]]
        A list of tuples that map ontology entries to one another
    symmetric : Optional[bool]
        If True, the mappings are interpreted as symmetric and will be applied
        in both directions
    """
    def __init__(self, statements, mappings=None, symmetric=True):
        self.statements = statements
        if mappings is None:
            self.mappings = _load_default_mappings()
        else:
            self.mappings = mappings
        self.symmetric = symmetric
        if self.symmetric:
            self._add_reverse_map()

    def map_statements(self):
        """Run the ontology mapping on the statements."""
        for stmt in self.statements:
            agents = stmt.agent_list()
            for agent in stmt.agent_list():
                if agent is None:
                    continue
                all_mappings = []
                for db_name, db_id in agent.db_refs.items():
                    if db_name == 'UN':
                        db_id = db_id[0][0]
                    mappings = self._map_id(db_name, db_id)
                    all_mappings += mappings
                for map_db_name, map_db_id in all_mappings:
                    if map_db_name in agent.db_refs:
                        continue
                    if map_db_name == 'UN':
                        agent.db_refs['UN'] = [(map_db_id, 1.0)]
                    else:
                        agent.db_refs[map_db_name] = map_db_id

    def _add_reverse_map(self):
        for m1, m2 in self.mappings:
            if (m2, m1) not in self.mappings:
                self.mappings.append((m2, m1))

    def _map_id(self, db_name, db_id):
        mappings = []
        # TODO: This lookup should be optimized using a dict
        for m1, m2 in self.mappings:
            if m1 == (db_name, db_id) or \
                m1 == (db_name, db_id.lower()):
                mappings.append(m2)
        return mappings

def _load_default_mappings():
    return [(('EIDOS', 'entities/x'), ('BBN', 'entities/y'))]




def _load_wm_map():
    path_here = os.path.dirname(os.path.abspath(__file__))
    ontomap_file = os.path.join(path_here, '../resources/wm_ontomap.tsv')
    mappings = {}
    with open(ontomap_file, 'r') as fh:
        for line in fh.readlines():
            s, se, t, te, score = line.split('\t')
            if s == 'eidos':
                s = 'UN'
            else:
                s = s.upper()
            if t == 'eidos':
                t = 'UN'
            else:
                t = t.upper()
            if s == 'SOFIA':
                parts = se.split('/')[1:]
                se = '/'.join(parts)
            if t == 'SOFIA':
                parts = te.split('/')[1:]
                te = '/'.join(parts)
            if s == 'BBN':
                se = 'events/%s' % se
            if t == 'BBN':
                te = 'events/%s' % te
            if (s, se) in mappings:
                if mappings[(s, se)][1] < score:
                    mappings[(s, se)] = ((t, te), score)
            else:
                mappings[(s, se)] = ((t, te), score)
    ontomap = []
    for s, ts in mappings.items():
        ontomap.append((s, ts[0]))
    return ontomap

try:
    wm_ontomap = _load_wm_map()
except Exception as e:
    wm_ontomap = []
