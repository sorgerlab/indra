import os
import logging
from functools import lru_cache


logger = logging.getLogger(__name__)


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
            self.mappings = []
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
                for map_db_name, map_db_id, score, orig_db_name in all_mappings:
                    if map_db_name in agent.db_refs:
                        continue
                    if self.scored:
                        # If the original one is a scored grounding,
                        # we take that score and multiply it with the mapping
                        # score. Otherwise we assume the original score is 1.
                        try:
                            orig_score = agent.db_refs[orig_db_name][0][1]
                        except Exception:
                            orig_score = 1.0
                        agent.db_refs[map_db_name] = \
                            [(map_db_id, score * orig_score)]
                    else:
                        if map_db_name in {'WM', 'UN'}:
                            agent.db_refs[map_db_name] = [(map_db_id, 1.0)]
                        else:
                            agent.db_refs[map_db_name] = map_db_id

    def _add_reverse_map(self):
        for m1, m2 in self.mappings:
            if (m2, m1) not in self.mappings:
                self.mappings.append((m2, m1))

    @lru_cache(maxsize=100000)
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
                mappings.append((m2[0], m2[1], score, db_name))
        return mappings


def _load_wm_map(exclude_auto=None):
    """Load an ontology map for world models.

    exclude_auto : None or list[tuple]
        A list of ontology mappings for which automated mappings should be
        excluded, e.g. [(HUME, UN)] would result in not using mappings
        from HUME to UN.
    """
    exclude_auto = [] if not exclude_auto else exclude_auto
    path_here = os.path.dirname(os.path.abspath(__file__))
    ontomap_file = os.path.join(path_here, '../resources/wm_ontomap.tsv')
    mappings = {}

    def map_entry(reader, entry):
        """Remap the readers and entries to match our internal standards."""
        if reader == 'WM':
            namespace = 'WM'
            entry_id = entry
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
            # Skip automated mappings when they should be excluded
            if (s, t) not in exclude_auto:
                # We first do the forward mapping
                if (s, se, t) in mappings:
                    if mappings[(s, se, t)][1] < score:
                        mappings[(s, se, t)] = ((t, te), score)
                else:
                    mappings[(s, se, t)] = ((t, te), score)
            # Then we add the reverse mapping
            if (t, s) not in exclude_auto:
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
