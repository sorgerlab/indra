from indra.pysb_assembler import PysbAssembler
from indra.trips import trips_api
from os.path import dirname, join
import indra.statements

test_small_file = join(dirname(__file__), 'test_small.xml')

def test_trips_processor_online():
    """Smoke test to see if imports and executes without error. Doesn't
    check for correctness of parse or of assembled model."""
    pa = PysbAssembler()
    tp = trips_api.process_text('BRAF phosphorylates MEK1 at Ser222.')
    pa.add_statements(tp.statements)

def test_trips_processor_offline():
    """Smoke test to see if imports and executes without error. Doesn't
    check for correctness of parse or of assembled model."""
    pa = PysbAssembler()
    tp = trips_api.process_xml(open(test_small_file).read())
    pa.add_statements(tp.statements)
    model = pa.make_model()

def test_bind():
    txt = 'The receptor tyrosine kinase EGFR binds the growth factor ligand EGF.'
    tp = trips_api.process_text(txt)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(has_hgnc_ref(st.members[0])
    assert(has_hgnc_ref(st.members[1])

def test_complex_bind():
    txt = 'The EGFR-EGF complex binds another EGFR-EGF complex.'
    tp = trips_api.process_text(txt)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(has_hgnc_ref(st.members[0])
    assert(has_hgnc_ref(st.members[1])
    assert(st.members[0].bound_to is not None)
    assert(st.members[1].bound_to is not None)

def is_complex(statement):
    isinstance(statement, indra.statements.Complex)

def has_hgnc_ref(agent):
    return 'hgnc' in agent.db_refs.keys()
