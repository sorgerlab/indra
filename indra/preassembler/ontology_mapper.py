import os
import rdflib


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
    scored : Optional[bool]
        If True, the mappings are assumed to be scored and the scores are
        propagated into the mapped groundings. If False, the scores don't
        need to be provided in the mappings and even if they are, they are
        ignored. Default: False
    """
    def __init__(self, statements, mappings=None, symmetric=True,
                 scored=False):
        self.statements = statements
        if mappings is None:
            self.mappings = _load_default_mappings()
        else:
            self.mappings = mappings
        self.symmetric = symmetric
        if self.symmetric:
            self._add_reverse_map()
        self.scored = scored

    def map_statements(self):
        """Run the ontology mapping on the statements."""
        for stmt in self.statements:
            for agent in stmt.agent_list():
                if agent is None:
                    continue
                all_mappings = []
                for db_name, db_id in agent.db_refs.items():
                    if isinstance(db_id, list):
                        db_id = db_id[0][0]
                    mappings = self._map_id(db_name, db_id)
                    all_mappings += mappings
                for map_db_name, map_db_id, score in all_mappings:
                    if map_db_name in agent.db_refs:
                        continue
                    if self.scored:
                        agent.db_refs[map_db_name] = [(map_db_id, score)]
                    else:
                        if map_db_name in ('UN', 'HUME'):
                            agent.db_refs[map_db_name] = [(map_db_id, 1.0)]
                        else:
                            agent.db_refs[map_db_name] = map_db_id

    def _add_reverse_map(self):
        for m1, m2 in self.mappings:
            if (m2, m1) not in self.mappings:
                self.mappings.append((m2, m1))

    def _map_id(self, db_name, db_id):
        mappings = []
        # TODO: This lookup should be optimized using a dict
        for mapping in self.mappings:
            if self.scored:
                m1, m2, score = mapping
            else:
                m1, m2 = mapping[:2]
                score = 1.0
            if m1 == (db_name, db_id) or \
                ((not isinstance(m1, list)) and
                 (m1 == (db_name, db_id.lower()))):
                mappings.append((m2[0], m2[1], score))
        return mappings


def _load_default_mappings():
    return [(('UN', 'entities/x'), ('HUME', 'entities/y'))]


def _load_wm_map():
    path_here = os.path.dirname(os.path.abspath(__file__))
    ontomap_file = os.path.join(path_here, '../resources/wm_ontomap.tsv')
    mappings = {}

    def make_hume_prefix_map():
        hume_ont = os.path.join(path_here, '../sources/hume/hume_ontology.rdf')
        graph = rdflib.Graph()
        graph.parse(os.path.abspath(hume_ont), format='nt')
        entry_map = {}
        for node in graph.all_nodes():
            entry = node.split('#')[1]
            # Handle "event" and other top-level entries
            if '/' not in entry:
                entry_map[entry] = None
                continue
            parts = entry.split('/')
            prefix, real_entry = parts[0], '/'.join(parts[1:])
            entry_map[real_entry] = prefix
        return entry_map

    hume_prefix_map = make_hume_prefix_map()

    def add_hume_prefix(hume_entry):
        """We need to do this because the HUME prefixes are missing"""
        prefix = hume_prefix_map[hume_entry]
        return '%s/%s' % (prefix, hume_entry)

    def map_entry(reader, entry):
        """Remap the reader and entry strings to match our internal standards."""
        if reader == 'eidos':
            namespace = 'UN'
            entry_id = entry
        elif reader == 'BBN':
            namespace = 'HUME'
            entry = entry.replace(' ', '_')
            entry_id = add_hume_prefix(entry)
        elif reader == 'sofia':
            namespace = 'SOFIA'
            # First chop off the Event/Entity prefix
            parts = entry.split('/')[1:]
            # Now we split each part by underscore and capitalize
            # each piece of each part
            parts = ['_'.join([p.capitalize() for p in part.split('_')])
                     for part in parts]
            # Finally we stick the entry back together separated by slashes
            entry_id = '/'.join(parts)
        else:
            return reader, entry
        return namespace, entry_id

    with open(ontomap_file, 'r') as fh:
        for line in fh.readlines():
            # Get each entry from the line
            s, se, t, te, score = line.strip().split('\t')
            score = float(score)
            # Map the entries to our internal naming standards
            s, se = map_entry(s, se)
            t, te = map_entry(t, te)
            # We first do the forward mapping
            if (s, se, t) in mappings:
                if mappings[(s, se, t)][1] < score:
                    mappings[(s, se, t)] = ((t, te), score)
            else:
                mappings[(s, se, t)] = ((t, te), score)
            # Then we add the reverse mapping
            if (t, te, s) in mappings:
                if mappings[(t, te, s)][1] < score:
                    mappings[(t, te, s)] = ((s, se), score)
            else:
                mappings[(t, te, s)] = ((s, se), score)
    ontomap = []
    for s, ts in mappings.items():
        ontomap.append(((s[0], s[1]), ts[0], ts[1]))
    return ontomap


try:
    wm_ontomap = _load_wm_map()
except Exception as e:
    wm_ontomap = []
