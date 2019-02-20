import logging
from indra.databases import hgnc_client, uniprot_client
from indra.statements import Agent, Phosphorylation, Evidence, BioContext, \
    RefContext

logger = logging.getLogger(__name__)


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

        # This will be the default name. If we get a gene name, it will
        # override this rawtext name.
        raw_text = entity_info['entityText']
        name = raw_text

        # Get the db refs.
        refs = {'TEXT': raw_text}
        for id_dict in entity_info['entityId']:
            if id_dict['source'] == 'Entrez':
                refs['EGID'] = id_dict['idString']
                hgnc_id = hgnc_client.get_hgnc_from_entrez(id_dict['idString'])
                if hgnc_id is not None:
                    # Check against what we may have already inferred from
                    # UniProt. If it disagrees with this, let it be. Inference
                    # from Entrez isn't as reliable.
                    if 'HGNC' in refs.keys():
                        if refs['HGNC'] != hgnc_id:
                            logger.info("HGNC id for Entrez id {EGID} did not "
                                        "match id inferred from UniProt id "
                                        "{UP}.".format(**refs))
                    else:
                        refs['HGNC'] = hgnc_id
            elif id_dict['source'] == 'UniProt':
                refs['UP'] = id_dict['idString']
                gene_name = uniprot_client.get_gene_name(id_dict['idString'])
                if gene_name is not None:
                    name = gene_name
                    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                    if hgnc_id is not None:
                        # Check to see if we have a conflict with an HGNC id
                        # found from the Entrez id. If so, overwrite with this
                        # one, in which we have greater faith.
                        if 'HGNC' in refs.keys() and refs['HGNC'] != hgnc_id:
                            logger.info("HGNC id for Entrez id {EGID} did not "
                                        "match id inferred from UniProt id "
                                        "{UP}.".format(**refs))
                        refs['HGNC'] = hgnc_id
            elif id_dict['source'] == 'Tax':
                refs['TAX'] = id_dict['idString']
            else:
                logger.warning("Unhandled id type: {source}={idString}"
                               .format(**id_dict))

        raw_coords = (entity_info['charStart'], entity_info['charEnd'])
        return Agent(name, db_refs=refs), raw_coords

    def _get_site(self, site_id):
        if site_id is None:
            return None, None, None
        site_info = self._entity_dict[site_id]
        site_text = site_info['attribute'][0]['value']
        residue = site_text[0]
        site_parts = site_text.split('-')
        position = None
        if len(site_parts) == 2:
            position = site_parts[1]
        coords = (site_info['charStart'], site_info['charEnd'])
        return residue, position, coords

    def _get_evidence(self, trigger_id, args, agent_coords, site_coords):
        """Get the evidence using the info in the trigger entity."""
        trigger_info = self._entity_dict[trigger_id]

        # Get the sentence index from the trigger word.
        s_idx_set = {self._entity_dict[eid]['sentenceIndex']
                     for eid in args.values()}
        i_min = min(s_idx_set)
        i_max = max(s_idx_set)

        text = '. '.join(self._sentences[i_min:(i_max+1)]) + '.'

        s_start = self._sentence_starts[i_min]
        annotations = {
            'agents': {'coords': [_fix_coords(coords, s_start)
                                  for coords in agent_coords]},
            'trigger': {'coords': _fix_coords([trigger_info['charStart'],
                                               trigger_info['charEnd']],
                                              s_start)}
            }
        if site_coords:
            annotations['site'] = {'coords': _fix_coords(site_coords, s_start)}

        return Evidence(text_refs=self._text_refs.copy(), text=text,
                        source_api='rlimsp', pmid=self._text_refs['pmid'],
                        annotations=annotations)

    def get_statements(self):
        stmts = []
        for rel_key, rel_info in self._relations.items():
            # Turn the arguments into a dict.
            args = {e['role']: e['entity_duid'] for e in rel_info['argument']}
            entity_args = args.copy()

            # Remove some special cases.
            trigger_id = entity_args.pop('TRIGGER')
            site_id = entity_args.pop('SITE', None)

            # Get the entity ids.
            entities = {role: self._get_agent(eid)
                        for role, eid in entity_args.items()}

            rel_type = rel_info['relationType']
            if rel_type == 'PHOSPHORYLATION':

                # Get the agents.
                enz, enz_coords = entities.get('KINASE', (None, None))
                sub, sub_coords = entities.get('SUBSTRATE', (None, None))

                # Get the site
                residue, position, site_coords = self._get_site(site_id)

                # Get the evidence
                ev = self._get_evidence(trigger_id, args,
                                        [enz_coords, sub_coords],
                                        site_coords)

                # Turn taxonomy into context, sub TAX takes precedence
                tax = None
                if enz and 'TAX' in enz.db_refs:
                    tax = enz.db_refs.pop('TAX')
                if sub and 'TAX' in sub.db_refs:
                    tax = sub.db_refs.pop('TAX')
                if tax is not None:
                    context = \
                        BioContext(species=RefContext(tax,
                                                      {'TAXONOMY': tax}))
                    ev.context = context

                stmts.append(Phosphorylation(enz, sub, residue=residue,
                                             position=position,
                                             evidence=[ev]))
            else:
                logger.warning("Unhandled statement type: %s" % rel_type)

        return stmts


def _fix_coords(coords, offset):
    """Adjust the entity coordinates to the beginning of the sentence."""
    if coords is None:
        return None
    return tuple([n - offset for n in coords])
