import os
from indra.sources import eidos
from indra.statements import Influence


path_this = os.path.dirname(os.path.abspath(__file__))
test_json = os.path.join(path_this, 'eidos_test.json')


def test_process_json():
    ep = eidos.process_json_file(test_json)
    assert ep is not None
    assert len(ep.statements) == 1
    stmt = ep.statements[0]
    assert isinstance(stmt, Influence)
    assert stmt.subj_delta.get('polarity') == 1
    assert stmt.obj_delta.get('polarity') == -1
    assert stmt.subj_delta.get('adjectives') == ['large']
    assert stmt.obj_delta.get('adjectives') == ['seriously']
    print(stmt)
