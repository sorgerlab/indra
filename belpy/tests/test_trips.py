from belpy.pysb_assembler import PysbAssembler
from belpy.trips import trips_api


def test_trips_processor_online():
    """Smoke test to see if imports and executes without error. Doesn't
    check for correctness of parse or of assembled model."""
    pa = PysbAssembler()
    tp = trips_api.process_text('BRAF phosphorylates MEK1 at Ser222')
    pa.add_statements(tp.statements)

def test_trips_processor_offline():
    """Smoke test to see if imports and executes without error. Doesn't
    check for correctness of parse or of assembled model."""
    pa = PysbAssembler()
    tp = trips_api.process_xml(open('test_small.xml').read())
    pa.add_statements(tp.statements)
    model = pa.make_model()
