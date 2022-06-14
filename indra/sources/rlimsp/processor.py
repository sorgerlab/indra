import logging
import tqdm
from collections import Counter
from indra.statements.validate import assert_valid_statements
from indra.databases import hgnc_client, uniprot_client
from indra.statements import Agent, Phosphorylation, Autophosphorylation, \
    Evidence, BioContext, RefContext, get_valid_residue, \
    InvalidResidueError, MutCondition

logger = logging.getLogger(__name__)


class RlimspProcessor(object):
    """Convert RLIMS-P JSON into INDRA Statements."""

    def __init__(self, rlimsp_json, doc_id_type=None):
        self._json = rlimsp_json
        self.statements = []
        self.doc_id_type = doc_id_type
        self.processed_texts = []
        return

    def extract_statements(self):
        """Extract the statements from the json."""
        for p_info in tqdm.tqdm(self._json, desc='Processing RLIMS-P JSON'):
            para = RlimspParagraph(p_info, self.doc_id_type)
            if para._text not in self.processed_texts:
                self.processed_texts.append(para._text)
                stmts = para.get_statements()
                assert_valid_statements(stmts)
                self.statements.extend(stmts)
        return


class RlimspParagraph(object):
    """An object that represents a single RLIMS-P Paragraph."""
    def __init__(self, p_info, doc_id_type):
        self._text = p_info['text']
        self._sentences = []
        self._sentence_starts = []
        for s in p_info['sentence']:
            start = s['charStart']
            stop = s['charEnd']
            self._sentences.append(self._text[start:stop])
            self._sentence_starts.append(start)
        if 'pmid' in p_info and 'pmcid' in p_info:
            self._text_refs = {n.upper(): p_info[n] for n in ['pmid', 'pmcid']
                               if p_info[n]}
        elif doc_id_type:
            self._text_refs = {doc_id_type.upper(): p_info['docId']}
        else:
            logger.info('Could not establish text refs for paragraph.')
            self._text_refs = {}
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
        return get_agent_from_entity_info(entity_info)

    def _get_site(self, site_id):
        def get_aa_code(residue_str):
            try:
                res = get_valid_residue(residue_str)
                return res
            except InvalidResidueError as e:
                logger.info('%s' % e)
                return None

        if site_id is None:
            return None, None, None
        site_info = self._entity_dict[site_id]
        site_text = site_info['attribute'][0]['value']
        site_parts = site_text.split('-')
        position = None
        if len(site_parts) == 2:
            residue_str = site_parts[0]
            residue = get_aa_code(residue_str)
            position = site_parts[1]
        else:
            residue = get_aa_code(site_text)
        coords = (site_info['charStart'], site_info['charEnd'])
        return residue, position, coords

    def _get_evidence(self, trigger_info, args, agent_coords, site_coords):
        """Get the evidence using the info in the trigger entity."""
        # Get the sentence index from the trigger word.
        s_idx_set = {self._entity_dict[eid]['sentenceIndex']
                     for eid in args.values()
                     if 'sentenceIndex' in self._entity_dict[eid]}
        if s_idx_set:
            i_min = min(s_idx_set)
            i_max = max(s_idx_set)

            text = '. '.join(self._sentences[i_min:(i_max+1)]) + '.'

            s_start = self._sentence_starts[i_min]
            annotations = {
                'agents': {'coords': [_fix_coords(coords, s_start)
                                      for coords in agent_coords]},
                'trigger': {'coords': _fix_coords([trigger_info['charStart'],
                                                   trigger_info['charEnd']],
                                                  s_start),
                            'text': trigger_info['entityText']}
                }
        else:
            logger.info('Unable to get sentence index')
            annotations = {}
            text = None
        if site_coords:
            annotations['site'] = {'coords': _fix_coords(site_coords, s_start)}

        return Evidence(text_refs=self._text_refs.copy(), text=text,
                        source_api='rlimsp', pmid=self._text_refs.get('PMID'),
                        annotations=annotations)

    def get_statements(self):
        stmts = []
        for rel_key, rel_info in self._relations.items():
            # Turn the arguments into a dict.
            args = {e['role']: e['entity_duid'] for e in rel_info['argument']}
            entity_args = args.copy()

            # Remove some special cases.
            trigger_id = entity_args.pop('TRIGGER')
            trigger_info = self._entity_dict[trigger_id]
            site_id = entity_args.pop('SITE', None)

            # Get the entity ids.
            entities = {role: self._get_agent(eid)
                        for role, eid in entity_args.items()}

            rel_type = rel_info['relationType']
            if rel_type == 'PHOSPHORYLATION':

                # Get the agents.
                enz, enz_coords = entities.get('KINASE', (None, None))
                sub, sub_coords = entities.get('SUBSTRATE', (None, None))
                if sub is None:
                    continue

                trigger_text = trigger_info.get('entityText')
                if enz is not None and enz.name == sub.name and \
                        'auto' in trigger_text:
                    is_autophos = True
                else:
                    is_autophos = False

                # Get the site
                residue, position, site_coords = self._get_site(site_id)

                # Get the evidence
                ev = self._get_evidence(trigger_info, args,
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

                if is_autophos:
                    stmt = Autophosphorylation(sub, residue=residue,
                                               position=position,
                                               evidence=[ev])
                else:
                    stmt = Phosphorylation(enz, sub, residue=residue,
                                           position=position,
                                           evidence=[ev])
                stmts.append(stmt)
            else:
                logger.warning("Unhandled statement type: %s" % rel_type)

        return stmts


def get_agent_from_entity_info(entity_info):
    """Return an INDRA Agent by processing an entity_info dict."""
    # This will be the default name. If we get a gene name, it will
    # override this rawtext name.
    raw_text = entity_info['entityText']
    name = raw_text

    # Get the db refs.
    refs = {'TEXT': raw_text}
    entries = entity_info['entityId']
    if entries is None or entries == {'$undefined': True}:
        entries = []
    ref_counts = Counter([entry['source'] for entry in entries])
    for source, count in ref_counts.items():
        if source in ('Entrez', 'UniProt') and count > 1:
            logger.info('%s has %d entries for %s, skipping'
                        % (raw_text, count, source))
            return None, None
    muts = []
    for id_dict in entries:
        if id_dict['source'] == 'Entrez':
            refs['EGID'] = id_dict['idString']
            hgnc_id = hgnc_client.get_hgnc_from_entrez(id_dict['idString'])
            if hgnc_id is not None:
                # Check against what we may have already inferred from
                # UniProt. If it disagrees with this, let it be. Inference
                # from Entrez isn't as reliable.
                if 'HGNC' in refs.keys():
                    if refs['HGNC'] != hgnc_id:
                        msg = ('HGNC:%s previously set does not'
                               ' match HGNC:%s from EGID:%s') % \
                               (refs['HGNC'], hgnc_id, refs['EGID'])
                        logger.info(msg)
                else:
                    refs['HGNC'] = hgnc_id
        elif id_dict['source'] == 'UniProt':
            refs['UP'] = id_dict['idString']
            hgnc_id = uniprot_client.get_hgnc_id(id_dict['idString'])
            if hgnc_id:
                # Check to see if we have a conflict with an HGNC id
                # found from the Entrez id. If so, overwrite with this
                # one, in which we have greater faith.
                if 'HGNC' in refs.keys() and refs['HGNC'] != hgnc_id:
                    msg = ('Inferred HGNC:%s from UP:%s does not'
                           ' match HGNC:%s from EGID:%s') % \
                          (refs['HGNC'], refs['UP'], hgnc_id,
                           refs['EGID'])
                    logger.info(msg)
                refs['HGNC'] = hgnc_id
                name = hgnc_client.get_hgnc_name(hgnc_id)
            else:
                gene_name = uniprot_client.get_gene_name(id_dict['idString'])
                if gene_name is not None:
                    name = gene_name
        elif id_dict['source'] in ('Tax', 'NCBI'):
            # Note that TAX is non-standard but it's popped out later in the
            # extraction process
            refs['TAX'] = id_dict['idString']
        elif id_dict['source'] == 'CHEBI':
            refs['CHEBI'] = 'CHEBI:%s' % id_dict['idString']
        # These we take as is
        elif id_dict['source'] in ('MESH', 'OMIM'):
            refs[id_dict['source']] = id_dict['idString']
        # CTD is sometimes used for MESH chemical IDs but can also be just '-'
        elif id_dict['source'] == 'CTD':
            if id_dict['idString'] != '-':
                refs['MESH'] = id_dict['idString']
        # Handle mutations
        elif id_dict['source'] == 'Unk' and \
                id_dict['entityType'] == 'ProteinMutation':
            # {'idString': 'p|SUB|Y|268|A', 'source': 'Unk',
            #  'tool': 'PubTator', 'entityType': 'ProteinMutation'}
            # Mpk1(Y268A)'
            if id_dict['idString'].startswith('p|SUB|'):
                try:
                    # Handle special cases like p|SUB|A|30|P;RS#:104893878
                    parts = id_dict['idString'].split(';')[0].split('|')
                    residue_from, pos, residue_to = parts[2:5]
                    mut = MutCondition(pos, residue_from, residue_to)
                    muts.append(mut)
                except Exception as e:
                    logger.info('Could not process mutation %s' %
                                id_dict['idString'])
            else:
                logger.info('Unhandled mutation: %s' % id_dict['idString'])
        else:
            logger.warning("Unhandled id type: {source}={idString}"
                           .format(**id_dict))

    raw_coords = (entity_info['charStart'], entity_info['charEnd'])
    return Agent(name, db_refs=refs, mutations=muts), raw_coords


def _fix_coords(coords, offset):
    """Adjust the entity coordinates to the beginning of the sentence."""
    if coords is None:
        return None
    return tuple([n - offset for n in coords])
