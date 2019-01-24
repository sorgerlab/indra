from indra.statements import Complex
from indra.sources.hprd import HprdProcessor
from indra.sources.hprd.api import HprdComplexRow

test_complex_rows = [
        HprdComplexRow(HPRD_COMPLEX_ID='COM_1', HPRD_PROTEIN_ID='00011',
                       HGNC_SYMBOL='ASCL1', REFSEQ_ID='NP_004307.2',
                       EXPT_TYPES='in vitro;in vivo', PMIDS='8900141,8948587'),
        HprdComplexRow(HPRD_COMPLEX_ID='COM_1', HPRD_PROTEIN_ID='00918',
                       HGNC_SYMBOL='TCF3', REFSEQ_ID='NP_003191.1',
                       EXPT_TYPES='in vitro;in vivo', PMIDS='8900141,8948587'),
        HprdComplexRow(HPRD_COMPLEX_ID='COM_1', HPRD_PROTEIN_ID='02809',
                       HGNC_SYMBOL='MEF2C', REFSEQ_ID='NP_002388.2',
                       EXPT_TYPES='in vitro;in vivo', PMIDS='8900141,8948587')]

def test_process_complexes():
    hp = HprdProcessor(test_complex_rows)
    assert isinstance(hp.statements, list)
    assert len(hp.statements) == 1
    stmt = hp.statements[0]
    assert isinstance(stmt, Complex)
    assert len(stmt.members) == 3
    assert set([ag.name for ag in stmt.members]) == {'ASCL1', 'TCF3', 'MEF2C'}

if __name__ == '__main__':
    test_process_complexes()
