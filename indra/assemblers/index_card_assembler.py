import json
from indra.statements import *
from indra.literature import id_lookup

global_submitter = 'sorgerlab'

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

def assemble_modification(stmt):
    card = IndexCard()
    card.card['pmc_id'] = get_pmc_id(stmt)
    card.card['submitter'] = global_submitter
    card.card['evidence'] = get_evidence_text(stmt)
    card.card['interaction']['interaction_type'] = 'adds_modification'
    card.card['interaction']['modifications'] = [{
            'feature_type': 'modification_feature',
            'modification_type': stmt.__class__.__name__.lower(),
            'location': stmt.position,
            'aa_code': stmt.residue
            }
        ]
    return card

def get_pmc_id(stmt):
    for ev in stmt.evidence:
        pmc_id = id_lookup(ev.pmid)['pmcid']
    if pmc_id.startswith('PMC'):
        pmc_id = pmc_id[3:]
    return pmc_id

def get_evidence_text(stmt):
    ev_txts = [ev.text for ev in stmt.evidence]
    return ev_txts
