import logging

from indra.statements import Agent, Phosphorylation, Evidence

logger = logging.getLogger("rlimsp_processor")


class RlimspProcessor(object):
    """Convert RLIMS-P JSON into INDRA Statements."""

    def __init__(self, rlimsp_json):
        self._json = rlimsp_json
        self.entity_dict = {}
        self.statements = []
        self.extract_statements()
        return

    def _get_agent(self, entity_id):
        """Convert the entity dictionary into an INDRA Agent."""
        info = self.entity_dict.get(entity_id)
        if info is None:
            logger.warning("Entity key did not resolve to entity.")
            return None

        # Keep it simple for now.
        name = info['entityText']

        # Get the db refs.
        refs = {'TEXT': name}
        for id_dict in info['entityId']:
            if id_dict['source'] == 'Entrez':
                refs['EGID'] = id_dict['idString']
            elif id_dict['source'] == 'UniProt':
                refs['UP'] = id_dict['idString']
            else:
                # FIXME: for development purposes only.
                print("Unhandled id type: {source}={idString}"
                      .format(**id_dict))

        return Agent(name, db_refs=refs)

    def _make_statement(self, rel_dict, text, text_refs):
        """Extract the information necessary to make an INDRA Statement."""
        stmt = None
        rel_type = rel_dict['relationType']
        entities = {e['role']: e['entity_duid'] for e in rel_dict['argument']}
        ev = Evidence(text=text, source_api='RLIMS-P',
                      text_refs=text_refs)
        if rel_type == 'PHOSPHORYLATION':
            enz = self._get_agent(entities.get('KINASE'))
            sub = self._get_agent(entities.get('SUBSTRATE'))
            stmt = Phosphorylation(enz, sub, evidence=[ev])
        else:
            # FIXME: For dev purposes only
            print("Unhandled statement type: %s" % rel_type)
        return stmt

    def extract_statements(self):
        """Extract the statements from the json."""
        for extraction in self._json:
            self.entity_dict = extraction['entity']
            text = extraction['text']
            text_refs = {'pmid': extraction['pmid'],
                         'pmcid': extraction['pmcid']}
            for rel_key, rel_info in extraction['relation'].items():
                stmt = self._make_statement(rel_info, text, text_refs.copy())
                self.statements.append(stmt)
        return
