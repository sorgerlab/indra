import pandas as pd
from collections import namedtuple
from indra.sources.hprd.processor import HprdProcessor

_hprd_id_cols = ['HPRD_ID', 'HGNC_SYMBOL', 'REFSEQ_GENE', 'REFSEQ_PROTEIN',
                 'EGID', 'OMIM', 'UNIPROT', 'NAME']

_cplx_cols = ['CPLX_ID', 'HPRD_ID', 'HGNC_SYMBOL', 'REFSEQ_PROTEIN',
              'EVIDENCE', 'PMIDS']

_ptm_cols = ['HPRD_ID', 'HGNC_SYMBOL', 'HPRD_ISOFORM', 'REFSEQ_PROTEIN',
             'POSITION', 'RESIDUE', 'ENZ_HGNC_SYMBOL', 'ENZ_HPRD_ID',
             'MOD_TYPE', 'EVIDENCE',
             'PMIDS']


def process_from_flat_files(id_mappings_file, complexes_file=None,
                            ptm_file=None):
    id_df = pd.read_csv(id_mappings_file, delimiter='\t', names=_hprd_id_cols,
                        dtype='str')
    id_df = id_df.set_index('HPRD_ID')
    if complexes_file is None and  ptm_file is None:
        raise ValueError('At least one of complexes_file or '
                         'ptm_file must be specified.')
    cplx_df = None
    if complexes_file:
        cplx_df = pd.read_csv(complexes_file, delimiter='\t', names=_cplx_cols,
                              dtype='str')
    ptm_df = None
    if ptm_file:
        ptm_df = pd.read_csv(ptm_file, delimiter='\t', names=_ptm_cols,
                             dtype='str', na_values='-')
    hprd_processor = HprdProcessor(id_df, cplx_df, ptm_df)
    return hprd_processor

