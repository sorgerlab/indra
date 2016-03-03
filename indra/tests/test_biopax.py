from indra.java_vm import autoclass, cast
from indra.biopax import biopax_api
import indra.biopax.processor as bpc
from indra.pysb_assembler import PysbAssembler

bp = biopax_api.process_owl('data/biopax_test.owl')
uri_prefix = 'http://purl.org/pc2/7/'

'''
def test_hyphenated_agent_names():
    """This query should contain reactions with agent names RAF1-BRAF,
    which need to be canonicalized to Python-compatible names before
    model assembly."""
    bp.get_phosphorylation()
    pa = PysbAssembler()
    pa.add_statements(bp.statements)
    pa.make_model()
'''

def test_paxtools_autoclass():
    autoclass('org.biopax.paxtools.impl.level3.ProteinImpl')

def test_biopaxpattern_autoclass():
    autoclass('org.biopax.paxtools.pattern.PatternBox')

def test_cpath_autoclass():
    autoclass('cpath.client.CPathClient')

def test_listify():
    assert(bpc.listify(1) == [1])
    assert(bpc.listify([1,2] == [1,2]))
    assert(bpc.listify([1] == [1]))

def test_list_listify():
    assert(bpc.list_listify([1]) == [[1]])
    assert(bpc.list_listify([1,2]) == [[1],[2]])
    assert(bpc.list_listify([1, [1,2]]) == [[1], [1,2]])

def test_get_combinations():
    combs = [c for c in bpc.get_combinations([1, 2])]
    assert(combs == [(1,2)])
    combs = [c for c in bpc.get_combinations([1, [3,4]])]
    assert(combs == [(1,3), (1,4)])

def test_has_members_er():
    bpe = bp.model.getByID(uri_prefix +\
                     'ProteinReference_971cec47bcd850e2b7d602f0416edacf')
    bpe = cast(bpc.bp('ProteinReference'), bpe)
    assert(bpc.has_members(bpe))
    
    bpe = bp.model.getByID('http://identifiers.org/uniprot/P56159')
    bpe = cast(bpc.bp('ProteinReference'), bpe)
    assert(not bpc.has_members(bpe))

def test_has_members_pe():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_117345.2')
    bpe = cast(bpc.bp('Protein'), bpe)
    assert(bpc.has_members(bpe))

def test_has_members_pe2():
    bpe = bp.model.getByID(uri_prefix + 'Protein_7d526475fd43d0a07ca1a596fe81aae0')
    bpe = cast(bpc.bp('Protein'), bpe)
    assert(not bpc.has_members(bpe))

def test_is_pe():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_117345.2')
    bpe = cast(bpc.bp('Protein'), bpe)
    assert(bpc.is_physical_entity(bpe))
    
def test_is_pe2():
    bpe = bp.model.getByID(uri_prefix +\
                     'ProteinReference_971cec47bcd850e2b7d602f0416edacf')
    bpe = cast(bpc.bp('ProteinReference'), bpe)
    assert(not bpc.is_physical_entity(bpe))
