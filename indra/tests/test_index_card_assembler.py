import json
import jsonschema
from indra.statements import *
from indra.assemblers.index_card_assembler import *

schema_path = os.path.dirname(os.path.abspath(__file__)) +\
             '/../resources/index_card_schema.json'
schema = json.load(open(schema_path, 'rt'))

braf = Agent('BRAF', db_refs={'UP': 'P15056'})
map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
ev = Evidence(source_api='reach', text='BRAF phosphorylates MAP2K1.', pmid='22833081')
stmt_phos = Phosphorylation(braf, map2k1, 'S', '222', evidence=ev)

brafmut = Agent('BRAF', db_refs={'UP': 'P15056'},
                mods=[ModCondition('phosphorylation', 'S', '596')],
                mutations=[MutCondition('600', 'V', 'E')],
                bound_conditions=[BoundCondition(Agent('BRAF'), True)])
stmt_phos2 = Phosphorylation(brafmut, map2k1, evidence=ev)
stmt_dephos = Dephosphorylation(brafmut, map2k1, evidence=ev)

def test_get_pmc_id():
    pmc_id = get_pmc_id(stmt_phos)
    assert(pmc_id == 'PMC4849135')

def test_get_evidence_text():
    ev_txt = get_evidence_text(stmt_phos)
    assert(len(ev_txt) == 1)
    assert(ev_txt[0] == 'BRAF phosphorylates MAP2K1.')

def test_assemble_phosphorylation():
    card = assemble_modification(stmt_phos)
    jsonschema.validate(card.card, schema)
    print card.get_string()
    print

def test_assemble_phosphorylation2():
    card = assemble_modification(stmt_phos2)
    jsonschema.validate(card.card, schema)
    print card.get_string()
    print

def test_assemble_dephosphorylation():
    card = assemble_modification(stmt_dephos)
    jsonschema.validate(card.card, schema)
    print card.get_string()
    print

def test_assemble_complex():
    card = assemble_complex(Complex([braf, brafmut, map2k1], evidence=ev))
    jsonschema.validate(card.card, schema)
    print card.get_string()
    print

def test_assemble_multiple():
    ia = IndexCardAssembler()
    ia.add_statements([stmt_phos, stmt_dephos])
    ia.make_model()
    ia.print_model()
    ia.save_model('/dev/null')

def test_get_participant():
    participant = get_participant(brafmut)
    print participant
