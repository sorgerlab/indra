import json
from indra.statements import *
from indra.literature import id_lookup
from indra.databases import hgnc_client, uniprot_client

global_submitter = 'cure'

class IndexCardAssembler(object):
    def __init__(self, statements=None):
        if statements is None:
            self.statements =  []
        else:
            self.statements = statements
        self.cards = []

    def add_statements(self, statements):
        self.statements.extend(statements)

    def make_model(self):
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                assemble_modification(stmt)
            elif isinstance(stmt, Complex):
                assemble_complex(stmt)
            else:
                print 'Assembly not defined for %s' % type(stmt)

class IndexCard(object):
    def __init__(self):
        self.card  = {
            'pmc_id': None,
            'submitter': None,
            'interaction': {
                'negative_information': False,
                'interaction_type': None,
                'participant_a': {
                    'entity_type': None,
                    'entity_text': None,
                    'identifier': None
                    },
                'participant_b': {
                    'entity_type': None,
                    'entity_text': None,
                    'identifier': None
                    }
                }
            }

    def get_string(self):
        return json.dumps(self.card)

def assemble_complex(stmt):
    card = IndexCard()
    card.card['pmc_id'] = get_pmc_id(stmt)
    card.card['submitter'] = global_submitter
    card.card['evidence'] = get_evidence_text(stmt)
    card.card['interaction']['interaction_type'] = 'complexes_with'
    card.pop('participant_b', None)
    # NOTE: fill out entity_text
    card.card['participant_a']['entity_type'] = 'complex'
    card.card['participant_a']['entities'] = []
    for m in stmt.members:
        p = get_participant(m)
        card.card['participant_a']['entities'].append(p)


def assemble_modification(stmt):
    card = IndexCard()
    card.card['pmc_id'] = get_pmc_id(stmt)
    card.card['submitter'] = global_submitter
    card.card['evidence'] = get_evidence_text(stmt)
    mod_type = stmt.__class__.__name__.lower()
    if mod_type.startswith('de'):
        card.card['interaction']['interaction_type'] = 'removes_modification'
    else:
        card.card['interaction']['interaction_type'] = 'adds_modification'
    card.card['interaction']['modifications'] = [{
            'feature_type': 'modification_feature',
            'modification_type': stmt.__class__.__name__.lower(),
            }]
    if stmt.position is not None:
        pos = int(stmt.position)
        card.card['interaction']['modifications'][0]['location'] = pos
    if stmt.residue is not None:
        card.card['interaction']['modifications'][0]['aa_code'] =  stmt.residue

    card.card['interaction']['participant_a'] = get_participant(stmt.enz)
    card.card['interaction']['participant_b'] = get_participant(stmt.sub)
    return card

def get_participant(agent):
    # Handle missing Agent as generic protein
    if agent is None:
        participant = {
            'entity_text': [''],
            'entity_type': 'protein',
            'identifier': 'GENERIC'
            }
        return participant
    # The Agent is not missing
    participant = {}
    participant['entity_text'] = [agent.name]
    hgnc_id = agent.db_refs.get('HGNC')
    uniprot_id = agent.db_refs.get('UP')
    chebi_id = agent.db_refs.get('CHEBI')
    # If HGNC grounding is available, that is the first choice
    if hgnc_id:
        uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    if uniprot_id:
        uniprot_mnemonic = uniprot_client.get_mnemonic(uniprot_id)
        participant['identifier'] = 'UNIPROT:%s' % uniprot_mnemonic
        participant['entity_type'] = 'protein'
    elif chebi_id:
        # NOTE: we need to convert to PubChem ID
        participant['identifier'] = 'CHEBI:%s' % chebi_id
        participant['entity_type'] = 'chemical'
    else:
        participant['identifier'] = None
        participant['entity_type'] = None

    features = []
    not_features = []
    # Binding features
    for bc in agent.bound_conditions:
        feature = {
            'feature_type': 'binding_feature',
            'bound_to': {
                # NOTE: get type and grounding for bound to protein
                'entity_type': 'protein',
                'entity_text': [bc.agent.name],
                'identifier': None
                }
            }
        if bc.is_bound:
            features.append(feature)
        else:
            not_features.append(feature)
    # Modification features
    for mc in agent.mods:
        feature = {
            'feature_type': 'modification_feature',
            'modification_type': mc.mod_type.lower(),
            }
        if mc.position is not None:
            pos = int(mc.position)
            feature['location'] = pos
        if mc.residue is not None:
            feature['aa_code'] = mc.residue
        if mc.is_modified:
            features.append(feature)
        else:
            not_features.append(feature)
    # Mutation features
    for mc in agent.mutations:
        feature = {
            'feature_type': 'mutation_feature',
            'from_aa': mc.residue_from,
            'to_aa': mc.residue_to
            }
        if mc.position is not None:
            pos = mc.position
            feature['location'] = pos
        features.append(feature)
    if features:
        participant['features'] = features
    if not_features:
        participant['not_features'] = not_features
    return participant

def get_pmc_id(stmt):
    for ev in stmt.evidence:
        pmc_id = id_lookup(ev.pmid)['pmcid']
    if not pmc_id.startswith('PMC'):
        pmc_id = 'PMC' + pmc_id
    return pmc_id

def get_evidence_text(stmt):
    ev_txts = [ev.text for ev in stmt.evidence]
    return ev_txts
