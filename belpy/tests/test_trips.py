from belpy.pysb_assembler import PysbAssembler
from belpy.trips import trips_api

def test_trips_processor():
    """Smoke test to see if imports and executes without error. Doesn't
    check for correctness of parse or of assembled model."""
    pa = PysbAssembler()
    tp = trips_api.process_text('BRAF phosphorylates MEK1 at Ser222')
    pa.add_statements(tp.belpy_stmts)
    model = pa.make_model()
