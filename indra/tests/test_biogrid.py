import os
from indra.statements import Complex
from indra.sources.biogrid import BiogridProcessor

this_dir = os.path.dirname(__file__)
test_file = os.path.join(this_dir, 'biogrid_tests_data/biogrid_test.txt')


def test_biogrid_tsv():
    # Download biogrid file form the web and process it
    bp = BiogridProcessor(test_file, physical_only=False)

    # There are 50 statements in that file
    statements = bp.statements
    assert len(statements) == 50, len(statements)

    # Any given statement should be a complex, with appropriate evidence
    # This one is a human protein complex
    s0 = statements[33]
    assert isinstance(s0, Complex)
    ev = s0.evidence[0]
    assert ev.source_api == 'biogrid'
    assert ev.text is None
    assert ev.pmid is not None

    # The first statement in the file involves MAP2K4 and FLNC
    assert s0.members[0].name == 'OTUB1', s0
    assert s0.members[1].name == 'DDX23', s0
    assert set(s0.members[0].db_refs) == {'HGNC', 'UP', 'EGID'}
    assert set(s0.members[1].db_refs) == {'HGNC', 'UP', 'EGID'}
