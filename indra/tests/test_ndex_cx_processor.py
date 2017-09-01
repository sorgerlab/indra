from indra.sources.ndex_cx import process_cx_file
from indra.sources.ndex_cx.processor import NdexCxProcessor

def test_process_cx_file():
    ncp = process_cx_file('merged_BRCA1_formatted.cx')
    assert isinstance(ncp, NdexCxProcessor)
