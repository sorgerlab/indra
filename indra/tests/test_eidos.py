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
    print(stmt)
