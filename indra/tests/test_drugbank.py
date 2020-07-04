import os
from indra.sources import drugbank

test_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'drugbank_sample.xml')


def test_drugbank_sample():
    dp = drugbank.process_xml(test_file)
    assert len(dp.statements) == 1
    stmt = dp.statements[0]
    assert len(stmt.evidence) == 6
    assert all(ev.pmid for ev in stmt.evidence)
    assert all(ev.source_api == 'drugbank' for ev in stmt.evidence)
    drug = stmt.subj
    assert drug.name == 'lepirudin'
    assert drug.db_refs['DRUGBANK'] == 'DB00001'
    assert drug.db_refs['CAS'] == '138068-37-8'

    target = stmt.obj
    assert target.name == 'F2'
    assert target.db_refs['HGNC'] == '3535'
    assert target.db_refs['UP'] == 'P00734'
    assert target.db_refs['DRUGBANKV4.TARGET'] == 'BE0000048'
