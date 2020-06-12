import os
from indra.statements import Complex
from indra.sources.biogrid import BiogridProcessor

this_dir = os.path.dirname(__file__)
test_file = os.path.join(this_dir, 'biogrid_tests_data/biogrid_test.txt')


def test_biogrid_tsv():
    # Download biogrid file form the web and process it
    bp = BiogridProcessor(test_file)

    # There are 50 statements in that file
    statements = bp.statements
    assert len(statements) == 50

    # Any given statement should be a complex, with appropriate evidence
    s0 = statements[0]
    assert isinstance(s0, Complex)
    ev = s0.evidence[0]
    assert ev.source_api == 'biogrid'
    assert ev.text is None
    assert ev.pmid is not None

    # The first statement in the file involves MAP2K4 and FLNC
    assert str(s0.members[0]) == 'MAP2K4()'
    assert str(s0.members[1]) == 'FLNC()'
