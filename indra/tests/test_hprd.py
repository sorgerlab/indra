from os.path import join, abspath, dirname
from indra.statements import Complex
from indra.sources import hprd


test_dir = join(abspath(dirname(__file__)), 'hprd_tests_data')

id_file = join(test_dir, 'HPRD_ID_MAPPINGS.txt')


def test_process_complexes():
    cplx_file = join(test_dir, 'PROTEIN_COMPLEXES.txt')
    hp = hprd.process_from_flat_files(id_file, complexes_file=cplx_file)
    assert isinstance(hp, hprd.HprdProcessor)
    assert isinstance(hp.statements, list)
    assert len(hp.statements) == 3
    s0 = hp.statements[0]
    assert isinstance(s0, Complex)
    assert len(s0.members) == 3
    assert set([ag.name for ag in s0.members]) == \
            set(['ASCL1', 'TCF3', 'MEF2C'])
    assert s0.members[0].db_refs == \
            {'HGNC': '738', 'UP': 'P50553', 'EGID': '429'}
    assert s0.members[1].db_refs == \
            {'HGNC': '11633', 'UP': 'P15923', 'EGID': '6929'}
    assert s0.members[2].db_refs == \
            {'HGNC': '6996', 'UP': 'Q06413', 'EGID': '4208'}
    assert len(s0.evidence) == 2
    assert s0.evidence[0].pmid == '8900141'
    assert s0.evidence[0].source_api == 'hprd'
    assert s0.evidence[0].source_id == ('http://hprd.org/interactions?'
                                     'hprd_id=00011&isoform_id=00011_1'
                                     '&isoform_name=')
    assert s0.evidence[1].pmid == '8948587'

if __name__ == '__main__':
    test_process_complexes()
