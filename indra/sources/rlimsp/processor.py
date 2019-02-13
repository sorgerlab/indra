import logging

from indra.statements import Agent, Phosphorylation, Evidence

logger = logging.getLogger("rlimsp_processor")


class RlimspProcessor(object):
    """Convert RLIMS-P JSON into INDRA Statements."""

    def __init__(self, rlimsp_json):
        self._json = rlimsp_json
        self.statements = []
        self.extract_statements()
        return

    def extract_statements(self):
        """Extract the statements from the json."""
        for p_info in self._json:
            para = RlimspParagraph(p_info)
            self.statements.extend(para.get_statements())
        return


class RlimspParagraph(object):
    """An object that represents a single RLIMS-P Paragraph."""
    def __init__(self, p_info):
        self._text = p_info['text']
        self._sentences = []
        self._sentence_starts = []
        for s in p_info['sentence']:
            start = s['charStart']
            stop = s['charEnd']
            self._sentences.append(self._text[start:stop])
            self._sentence_starts.append(start)
        self._text_refs = {n: p_info[n] for n in ['pmid', 'pmcid']}
        self._relations = p_info['relation']
        self._entity_dict = p_info['entity']
        return

    def _get_agent(self, entity_id):
        """Convert the entity dictionary into an INDRA Agent."""
        if entity_id is None:
            return None

        entity_info = self._entity_dict.get(entity_id)
        if entity_info is None:
            logger.warning("Entity key did not resolve to entity.")
            return None

        # Keep it simple for now.
        name = entity_info['entityText']

        # Get the db refs.
        refs = {'TEXT': name}
        for id_dict in entity_info['entityId']:
            if id_dict['source'] == 'Entrez':
                refs['EGID'] = id_dict['idString']
            elif id_dict['source'] == 'UniProt':
                refs['UP'] = id_dict['idString']
            else:
                # FIXME: for development purposes only.
                print("Unhandled id type: {source}={idString}"
                      .format(**id_dict))

        raw_coords = (entity_info['charStart'], entity_info['charEnd'])
        return Agent(name, db_refs=refs), raw_coords

    def get_evidence(self, args, agent_coords):
        """Get the evidence using the info in the trigger entity."""
        trigger_info = self._entity_dict[args['TRIGGER']]

        # Get the sentence index from the trigger word.
        sent_idx = trigger_info['sentenceIndex']

        # Check for other sentence indices.
        other_sent_idx = {self._entity_dict[eid]['sentenceIndex']
                          for eid in args.values()}
        other_sent_idx -= {sent_idx}
        if other_sent_idx:
            logger.warning("Found another sentence index among the "
                           "entities: %s." % other_sent_idx)

        s_start = self._sentence_starts[sent_idx]
        annotations = {
            'agents': {'coords': [_fix_coords(coords, s_start)
                                  for coords in agent_coords]},
            'trigger': {'coords': _fix_coords([trigger_info['charStart'],
                                               trigger_info['charEnd']],
                                              s_start)}
            }

        return Evidence(text_refs=self._text_refs.copy(),
                        text=self._sentences[sent_idx],
                        source_api='RLIMS-P', pmid=self._text_refs['pmid'],
                        annotations=annotations)

    def get_statements(self):
        stmts = []
        for rel_key, rel_info in self._relations.items():
            # Turn the arguments into a dict.
            args = {e['role']: e['entity_duid'] for e in rel_info['argument']}

            # Get the entity ids.
            entities = {role: self._get_agent(eid) for role, eid in args.items() if role != 'TRIGGER'}

            # Check to make sure we didn't loose any. Roles are presumed
            # unique, but that may not always be true.
            if len(entities) != len(rel_info['argument']):
                logger.warning("Lost some arguments: %s vs. %s"
                               % (entities, rel_info['arguments']))

            rel_type = rel_info['relationType']
            if rel_type == 'PHOSPHORYLATION':

                # Get the agents.
                enz, enz_coords = entities.get('KINASE', (None, None))
                sub, sub_coords = entities.get('SUBSTRATE', (None, None))

                # Get the evidence
                ev = self.get_evidence(args, [enz_coords, sub_coords])

                stmts.append(Phosphorylation(enz, sub, evidence=[ev]))
            else:
                # FIXME: For dev purposes only
                print("Unhandled statement type: %s" % rel_type)

        return stmts


def _fix_coords(coords, offset):
    """Adjust the entity coordinates to the beginning of the sentence."""
    return tuple([n - offset for n in coords])
