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


if __name__ == '__main__':
    hprd_path = '/Users/johnbachman/Dropbox/1johndata/Knowledge File/Biology/Research/Big Mechanism/protmapper_paper/data/hprd'

    from os.path import join
    pd.set_option('display.max_columns', 15)
    pd.set_option('expand_frame_repr', False)
    hp = process_from_flat_files(
                    join(hprd_path, 'HPRD_ID_MAPPINGS.txt'),
                    #complexes_file=join(hprd_path, 'PROTEIN_COMPLEXES.txt'),
        ptm_file=join(hprd_path, 'POST_TRANSLATIONAL_MODIFICATIONS.txt'))

