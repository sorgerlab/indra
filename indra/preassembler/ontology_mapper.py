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
                    if db_name == 'EIDOS':
                        db_id = db_id[0][0]
                    mappings = self._map_id(db_name, db_id)
                    all_mappings += mappings
                for map_db_name, map_db_id in all_mappings:
                    if map_db_name in agent.db_refs:
                        continue
                    if map_db_name == 'EIDOS':
                        agent.db_refs['EIDOS'] = [(map_db_id, 1.0)]
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
            if m1 == (db_name, db_id):
                mappings.append(m2)
        return mappings

def _load_default_mappings():
    return [(('EIDOS', 'entities/x'), ('BBN', 'entities/y'))]
