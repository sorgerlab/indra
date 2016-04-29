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
            'pmc_id': None
            'submitter': None
            'interaction': {
                'negative_information': False,
                'interaction_type': None,
                'participant_a': {
                    'entity_type': None,
                    'entity_text': None,
                    'identifier': None
                    }
                'participant_b': {
                    'entity_type': None,
                    'entity_text': None,
                    'identifier': None
                    }
                }
            }

def assemble_modification(stmt):
    card = IndexCard()
    card['pmc_id'] = get_pmc_id(stmt)
    card['submitter'] = global_submitter
    card['evidence'] = get_evidence_text(stmt)
    card['interaction']['interaction_type'] = 'adds_modification'
    card['interaction']['modifications'] = [
        {
        'feature_type': 'modification_feature',
        'modification_type': str(type(stmt)).lower(),
        'location': stmt.position,
        'aa_code': stmt.residue
        }
        ]

def get_pmc_id(stmt):
    for ev in stmt.evidence:
        pmc_id = id_lookup(ev.pmid)['pmc']
    return pmc_id

def get_evidence_text(stmt):
    ev_txts = [ev.text for ev in stmt.evidence]
    return ev_txts
