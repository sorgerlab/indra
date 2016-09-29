from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import jsonschema
from indra.statements import *
from indra.assemblers.index_card_assembler import *
from os.path import dirname, abspath, join

schema_path = join(dirname(abspath(__file__)),
                   '../resources/index_card_schema.json')
with open(schema_path, 'rt') as fh:
    schema = json.load(fh)

braf = Agent('BRAF', db_refs={'UP': 'P15056'})
map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
ev = Evidence(source_api='reach', text='BRAF phosphorylates MAP2K1.',
              pmid='22833081')
stmt_phos = Phosphorylation(braf, map2k1, 'S', '222', evidence=ev)

brafmut = Agent('BRAF', db_refs={'UP': 'P15056'},
                mods=[ModCondition('phosphorylation', 'S', '596')],
                mutations=[MutCondition('600', 'V', 'E')],
                bound_conditions=[BoundCondition(Agent('BRAF'), True)])
stmt_phos2 = Phosphorylation(brafmut, map2k1, evidence=ev)
stmt_dephos = Dephosphorylation(brafmut, map2k1, evidence=ev)
stmt_complex = Complex([Agent('HRAS', db_refs={'UP': 'P01112'}),
                        Agent('GTP', db_refs={'CHEBI': '57600'})],
                        evidence=ev)
ev2 = Evidence(source_api='reach', text='BRAF phosphorylates MAP2K1.',
               pmid='22833081', epistemics={'direct': False})
stmt_phos_indirect = Phosphorylation(brafmut, map2k1, evidence=ev2)
stmt_autophos = Autophosphorylation(brafmut, 'S', '564', evidence=ev)

def test_get_pmc_id():
    pmc_id = get_pmc_id(stmt_phos)
    assert(pmc_id == 'PMC4849135')

def test_get_evidence_text():
    ev_txt = get_evidence_text(stmt_phos)
    assert(len(ev_txt) == 1)
    assert(ev_txt[0] == 'BRAF phosphorylates MAP2K1.')

def test_assemble_phosphorylation():
    card = assemble_modification(stmt_phos)
    card.card['pmc_id'] = get_pmc_id(stmt_phos)
    print(card.get_string())
    print()
    jsonschema.validate(card.card, schema)

def test_assemble_phosphorylation2():
    card = assemble_modification(stmt_phos2)
    card.card['pmc_id'] = get_pmc_id(stmt_phos2)
    print(card.get_string())
    print()
    jsonschema.validate(card.card, schema)

def test_assemble_phosphorylation_indirect():
    card = assemble_modification(stmt_phos_indirect)
    card.card['pmc_id'] = get_pmc_id(stmt_phos_indirect)
    print(card.get_string())
    print()
    jsonschema.validate(card.card, schema)

def test_assemble_dephosphorylation():
    card = assemble_modification(stmt_dephos)
    card.card['pmc_id'] = get_pmc_id(stmt_dephos)
    print(card.get_string())
    print()
    jsonschema.validate(card.card, schema)

def test_assemble_autophosphorylation():
    card = assemble_selfmodification(stmt_autophos)
    card.card['pmc_id'] = get_pmc_id(stmt_autophos)
    print(card.get_string())
    print()
    jsonschema.validate(card.card, schema)

def test_assemble_complex():
    st_complex = Complex([braf, brafmut, map2k1], evidence=ev)
    card = assemble_complex(st_complex)
    card.card['pmc_id'] = get_pmc_id(st_complex)
    print(card.get_string())
    print()
    jsonschema.validate(card.card, schema)

def test_assemble_multiple():
    ia = IndexCardAssembler(pmc_override='PMC1234567')
    ia.add_statements([stmt_phos, stmt_dephos])
    ia.make_model()
    ia.print_model()
    ia.save_model('/dev/null')

def test_get_participant():
    participant = get_participant(brafmut)
    print(participant)

def test_chemical():
    card = assemble_complex(stmt_complex)
    card.card['pmc_id'] = get_pmc_id(stmt_complex)
    print(card.get_string())
    print()
    jsonschema.validate(card.card, schema)
