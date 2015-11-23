from indra.pysb_assembler import PysbAssembler
from indra.trips import trips_api
from os.path import dirname, join

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
