from __future__ import absolute_import, print_function, unicode_literals
import json
from builtins import dict, str
from indra.statements import *
from indra.assemblers.kami_assembler import KamiAssembler


mek = Agent('MAP2K1', db_refs={'HGNC': '6840'})
erk = Agent('MAPK1', db_refs={'UP': 'P28482'})
braf = Agent('BRAF', db_refs={'HGNC': '1097'})
mek_braf = Agent('MAP2K1', bound_conditions=[BoundCondition(braf, True)],
                 db_refs={'HGNC': '6840'})
mek_no_braf = Agent('MAP2K1', bound_conditions=[BoundCondition(braf, False)],
                 db_refs={'HGNC': '6840'})
mek_phos = Agent('MAP2K1', mods=[ModCondition('phosphorylation', None,
                                              None, True)],
                 db_refs={'HGNC': '6840'})
mek_phos2 = Agent('MAP2K1', mods=[ModCondition('phosphorylation', 'S',
                                               '222', True)],
                 db_refs={'HGNC': '6840'})
mek_phos3 = Agent('MAP2K1', mods=[ModCondition('phosphorylation', 'S',
                                               '222', False)],
                 db_refs={'HGNC': '6840'})

def test_complex_no_conditions():
    stmt = Complex([mek, erk])
    ka = KamiAssembler()
    ka.add_statements([stmt])
    model = ka.make_model()
    assert isinstance(model, dict)
    assert isinstance(model['graphs'], list)
    assert isinstance(model['typing'], list)
    graph_list = model['graphs']
    assert len(graph_list) == 3
    assert len(graph_list[1]['graph']['edges']) == 4
    assert len(graph_list[1]['graph']['nodes']) == 5
    print(json.dumps(model, indent=1))

def test_complex_bound_condition():
    stmt = Complex([mek_braf, erk])
    ka = KamiAssembler()
    ka.add_statements([stmt])
    model = ka.make_model()
    assert isinstance(model, dict)
    assert isinstance(model['graphs'], list)
    assert isinstance(model['typing'], list)
    graph_list = model['graphs']
    assert len(graph_list) == 3
    assert len(graph_list[1]['graph']['edges']) == 6
    assert len(graph_list[1]['graph']['nodes']) == 7
    print(json.dumps(model, indent=1))

def test_complex_not_bound_condition():
    stmt = Complex([mek_no_braf, erk])
    ka = KamiAssembler()
    ka.add_statements([stmt])
    model = ka.make_model()
    assert isinstance(model, dict)
    assert isinstance(model['graphs'], list)
    assert isinstance(model['typing'], list)
    graph_list = model['graphs']
    assert len(graph_list) == 3
    assert len(graph_list[1]['graph']['edges']) == 6
    assert len(graph_list[1]['graph']['nodes']) == 7
    print(json.dumps(model, indent=1))

def test_complex_mod_condition():
    meks = [mek_phos, mek_phos2, mek_phos3]
    for mek in meks:
        stmt = Complex([mek, erk])
        ka = KamiAssembler()
        ka.add_statements([stmt])
        model = ka.make_model()
        assert isinstance(model, dict)
        assert isinstance(model['graphs'], list)
        assert isinstance(model['typing'], list)
        graph_list = model['graphs']
        assert len(graph_list) == 3
        assert len(graph_list[1]['graph']['edges']) == 5
        assert len(graph_list[1]['graph']['nodes']) == 6
        print(json.dumps(model, indent=1))

def test_mod_one_step():
    stmt = Phosphorylation(mek, erk)
    ka = KamiAssembler()
    ka.add_statements([stmt])
    model = ka.make_model()
    print(json.dumps(model, indent=1))
    assert isinstance(model, dict)
    assert isinstance(model['graphs'], list)
    assert isinstance(model['typing'], list)
    graph_list = model['graphs']
    assert len(graph_list) == 2
    assert len(graph_list[1]['graph']['edges']) == 3
    assert len(graph_list[1]['graph']['nodes']) == 4

def test_demod_one_step():
    stmt = Deubiquitination(mek, erk, 'K', '123')
    ka = KamiAssembler()
    ka.add_statements([stmt])
    model = ka.make_model()
    print(json.dumps(model, indent=1))
    assert isinstance(model, dict)
    assert isinstance(model['graphs'], list)
    assert isinstance(model['typing'], list)
    graph_list = model['graphs']
    assert len(graph_list) == 2
    assert len(graph_list[1]['graph']['edges']) == 3
    assert len(graph_list[1]['graph']['nodes']) == 4

def test_unique_id():
    egf = Agent('EGF')
    egfr = Agent('EGFR', bound_conditions=BoundCondition(egf, True))
    stmt = Complex([egfr, egfr])
    ka = KamiAssembler()
    ka.add_statements([stmt])
    model = ka.make_model()
    print(json.dumps(model, indent=1))

