from collections import namedtuple

HprdComplexRow = namedtuple('HprdComplexRow',
                        ['HPRD_COMPLEX_ID', 'HPRD_PROTEIN_ID', 'HGNC_SYMBOL',
                         'REFSEQ_ID', 'EXPT_TYPES', 'PMIDS'])


def process_from_flat_files(complexes_file=None, ptm_file=None):
    if complexes_file is None and phosphorylation_file is None:
        raise ValueError('At least one of complexes_file or '
                         'phosphorylation_file must be specified.')
    hprd_rows = []
    if ptm_file:
        pass
    if complexes_file:
        pass
