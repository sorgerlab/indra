from indra.statements import *
from indra.assemblers.index_card_assembler import *

braf = Agent('BRAF', db_refs={'UP': 'P15056'})
map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
ev = Evidence(source_api='reach', text='BRAF phosphorylates MAP2K1.', pmid='22833081')
stmt_phos = Phosphorylation(braf, map2k1, 'S', '222', evidence=ev)

brafmut = Agent('BRAF', db_refs={'UP': 'P15056'},
                mods=[ModCondition('phosphorylation', 'S', '596')],
                mutations=[MutCondition('600', 'V', 'E')],
                bound_conditions=[BoundCondition(Agent('BRAF'), True)])

def test_get_pmc_id():
    pmc_id = get_pmc_id(stmt_phos)
    assert(pmc_id == '4849135')

def test_get_evidence_text():
    ev_txt = get_evidence_text(stmt_phos)
    assert(len(ev_txt) == 1)
    assert(ev_txt[0] == 'BRAF phosphorylates MAP2K1.')

def test_assemble_phosphorylation():
    card = assemble_modification(stmt_phos)
    print card.get_string()

def test_get_participant():
    participant = get_participant(brafmut)
    print participant
