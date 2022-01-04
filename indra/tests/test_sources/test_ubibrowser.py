import os
from indra.statements import Ubiquitination, Deubiquitination
from indra.sources import ubibrowser
from indra.statements.validate import assert_valid_statements


e3_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'resources', 'ubibrowser_e3.txt')
dub_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'resources', 'ubibrowser_dub.txt')


def test_extract_statements():
    up = ubibrowser.process_file(e3_file, dub_file)
    assert len(up.statements) == 9
    assert isinstance(up.statements[0], Ubiquitination)
    assert isinstance(up.statements[-1], Deubiquitination)

    assert_valid_statements(up.statements)

    #1	AMFR2_HUMAN	A1AT_HUMAN	Q9UKV5	P01009	AMFR	SERPINA1
    # MEDLINE	16979136	Here we report that ...	RING	3	E3	H.sapiens
    e3_stmt = up.statements[0]
    assert e3_stmt.enz.name == 'AMFR'
    assert e3_stmt.enz.db_refs['UP'] == 'Q9UKV5'
    assert e3_stmt.sub.name == 'SERPINA1'
    assert e3_stmt.sub.db_refs['UP'] == 'P01009'
    assert len(e3_stmt.evidence) == 1
    assert e3_stmt.evidence[0].source_api == 'ubibrowser'
    assert e3_stmt.evidence[0].pmid == '16979136'
    assert e3_stmt.evidence[0].text.startswith('Here we report that')

    # 677	UBP33_HUMAN	ARRB1_HUMAN	Q8TEY7	P49407	USP33	ARRB1
    # MEDLINE	19363159	We now report the discovery that...	USP	1	DUB	H.sapiens
    dub_stmt = up.statements[-1]
    assert dub_stmt.enz.name == 'USP33'
    assert dub_stmt.enz.db_refs['UP'] == 'Q8TEY7'
    assert dub_stmt.sub.name == 'ARRB1'
    assert dub_stmt.sub.db_refs['UP'] == 'P49407'
    assert len(dub_stmt.evidence) == 1
    assert dub_stmt.evidence[0].source_api == 'ubibrowser'
    assert dub_stmt.evidence[0].pmid == '19363159'
    assert dub_stmt.evidence[0].text.startswith('We now report the discovery that')
