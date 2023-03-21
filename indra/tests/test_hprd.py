from os.path import join, abspath, dirname
import pytest
from indra.statements import Complex, Phosphorylation
from indra.sources import hprd


test_dir = join(abspath(dirname(__file__)), 'hprd_tests_data')


id_file = join(test_dir, 'HPRD_ID_MAPPINGS.txt')


def test_process_complexes():
    cplx_file = join(test_dir, 'PROTEIN_COMPLEXES.txt')
    hp = hprd.process_flat_files(id_file, complexes_file=cplx_file)
    assert isinstance(hp, hprd.HprdProcessor)
    assert isinstance(hp.statements, list)
    assert len(hp.statements) == 3
    s0 = hp.statements[0]
    assert isinstance(s0, Complex)
    assert len(s0.members) == 3
    assert set([ag.name for ag in s0.members]) == \
            set(['ASCL1', 'TCF3', 'MEF2C'])
    assert s0.members[0].db_refs == \
            {'HGNC': '738', 'UP': 'P50553', 'EGID': '429',
             'REFSEQ_PROT': 'NP_004307.2'}
    assert s0.members[1].db_refs == \
            {'HGNC': '11633', 'UP': 'P15923', 'EGID': '6929',
             'REFSEQ_PROT': 'NP_003191.1'}
    assert s0.members[2].db_refs == \
            {'HGNC': '6996', 'UP': 'Q06413', 'EGID': '4208',
             'REFSEQ_PROT': 'NP_002388.2'}
    assert len(s0.evidence) == 2
    assert s0.evidence[0].pmid == '8900141'
    assert s0.evidence[0].source_api == 'hprd'
    assert s0.evidence[0].annotations['evidence'] == ['in vivo']
    assert s0.evidence[0].source_id == ('http://hprd.org/interactions?'
                                     'hprd_id=00011&isoform_id=00011_1'
                                     '&isoform_name=Isoform_1')
    assert s0.evidence[1].pmid == '8948587'


def test_process_ptms():
    ptm_file = join(test_dir, 'POST_TRANSLATIONAL_MODIFICATIONS.txt')
    seq_file = join(test_dir, 'PROTEIN_SEQUENCES.txt')
    hp = hprd.process_flat_files(id_file, ptm_file=ptm_file, seq_file=seq_file)
    assert isinstance(hp, hprd.HprdProcessor)
    assert isinstance(hp.statements, list)
    assert len(hp.statements) == 13
    s0 = hp.statements[0]
    assert isinstance(s0, Phosphorylation)
    assert s0.enz.name == 'MAPK1'
    assert s0.enz.db_refs == {'UP': 'P28482', 'HGNC': '6871', 'EGID': '5594',
                              'REFSEQ_PROT': 'NP_002736.3'}
    assert s0.sub.name == 'TCF3'
    assert s0.sub.db_refs == {'UP': 'P15923', 'HGNC': '11633', 'EGID': '6929',
                              'REFSEQ_PROT': 'NP_003191.1'}
    assert s0.residue == 'T'
    assert s0.position == '355'
    assert len(s0.evidence) == 1
    assert s0.evidence[0].pmid == '14592976'
    assert s0.evidence[0].source_api == 'hprd'
    assert s0.evidence[0].annotations['evidence'] == ['in vivo']
    assert s0.evidence[0].annotations['site_motif'] == \
                    {'motif': 'NFSSSPSTPVGSPQG', 'respos': 8,
                     'off_by_one': False}


def test_process_ppis():
    ppi_file = join(test_dir, 'BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt')
    hp = hprd.process_flat_files(id_file, ppi_file=ppi_file)
    assert isinstance(hp, hprd.HprdProcessor)
    assert isinstance(hp.statements, list)
    assert len(hp.statements) == 5
    s0 = hp.statements[0]
    assert isinstance(s0, Complex)
    assert len(s0.members) == 2
    assert set([ag.name for ag in s0.members]) == set(['ITGA7', 'CHRNA1'])
    assert s0.members[0].db_refs == \
            {'HGNC': '6143', 'UP': 'Q13683', 'EGID': '3679',
             'REFSEQ_PROT': 'NP_001138468.1'}
    assert s0.members[1].db_refs == \
            {'HGNC': '1955', 'UP': 'P02708', 'EGID': '1134',
             'REFSEQ_PROT': 'NP_001034612.1'}
    assert len(s0.evidence) == 1
    assert s0.evidence[0].pmid == '10910772'
    assert s0.evidence[0].source_api == 'hprd'
    assert s0.evidence[0].annotations['evidence'] == ['in vivo']
    assert s0.evidence[0].source_id == ('http://hprd.org/interactions?'
                                     'hprd_id=02761&isoform_id=02761_1'
                                     '&isoform_name=Isoform_1')


def test_process_ptms_no_seq():
    with pytest.raises(ValueError):
        ptm_file = join(test_dir, 'POST_TRANSLATIONAL_MODIFICATIONS.txt')
        hp = hprd.process_flat_files(id_file, ptm_file=ptm_file)

