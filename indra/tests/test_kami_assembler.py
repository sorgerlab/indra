from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers.kami_assembler import KamiAssembler


mek = Agent('MAP2K1', db_refs={'HGNC': '6840'})
erk = Agent('MAPK1', db_refs={'UP': 'P28482'})

def test_phosphorylation_no_site():
    stmt = Complex([mek, erk])
    ka = KamiAssembler()
    ka.add_statements([stmt])
    model = ka.make_model()
    assert model == {'foo': 'bar'}

if __name__ == '__main__':
    test_phosphorylation_no_site()
