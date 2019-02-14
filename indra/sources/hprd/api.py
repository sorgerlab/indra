import pandas as pd
from collections import namedtuple
from indra.sources.hprd.processor import HprdProcessor

_hprd_id_cols = ['HPRD_ID', 'HGNC_SYMBOL', 'REFSEQ_GENE', 'REFSEQ_PROTEIN',
                 'EGID', 'OMIM', 'UNIPROT', 'NAME']

_cplx_cols = ['CPLX_ID', 'HPRD_ID', 'HGNC_SYMBOL', 'REFSEQ_PROTEIN',
              'EVIDENCE', 'PMIDS']

"""
_cplx_dtype = {
    'CPLX_ID': 'str',
    'HPRD_ID': 'str',
    'HGNC_SYMBOL': 'str',
    'REFSEQ_PROTEIN': 'str',
    'EVIDENCE': 'str',
    'PMIDS': 'str',
    }
"""

def process_from_flat_files(id_mappings_file, complexes_file=None,
                            ptm_file=None):
    id_df = pd.read_csv(id_mappings_file, delimiter='\t', names=_hprd_id_cols,
                        dtype='str')
    id_df = id_df.set_index('HPRD_ID')
    if complexes_file is None and phosphorylation_file is None:
        raise ValueError('At least one of complexes_file or '
                         'phosphorylation_file must be specified.')
    if complexes_file:
        cplx_df = pd.read_csv(complexes_file, delimiter='\t', names=_cplx_cols,
                              dtype='str')
    hprd_processor = HprdProcessor(id_df, cplx_df)
    return hprd_processor

if __name__ == '__main__':
    from os.path import join
    pd.set_option('display.max_columns', 15)
    pd.set_option('expand_frame_repr', False)
    hp = process_from_flat_files(
                            join(hprd_path, 'HPRD_ID_MAPPINGS.txt'),
                            join(hprd_path, 'PROTEIN_COMPLEXES.txt'))

