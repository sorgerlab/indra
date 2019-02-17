import pandas as pd
from collections import namedtuple
from protmapper.uniprot_client import load_fasta_sequences
from indra.sources.hprd.processor import HprdProcessor

_hprd_id_cols = ['HPRD_ID', 'HGNC_SYMBOL', 'REFSEQ_GENE', 'REFSEQ_PROTEIN',
                 'EGID', 'OMIM', 'UNIPROT', 'NAME']

_cplx_cols = ['CPLX_ID', 'HPRD_ID', 'HGNC_SYMBOL', 'REFSEQ_PROTEIN',
              'EVIDENCE', 'PMIDS']

_ptm_cols = ['HPRD_ID', 'HGNC_SYMBOL', 'HPRD_ISOFORM', 'REFSEQ_PROTEIN',
             'POSITION', 'RESIDUE', 'ENZ_HGNC_SYMBOL', 'ENZ_HPRD_ID',
             'MOD_TYPE', 'EVIDENCE',
             'PMIDS']

_ppi_cols = ['HGNC_SYMBOL_A', 'HPRD_ID_A', 'REFSEQ_PROTEIN_A',
             'HGNC_SYMBOL_B', 'HPRD_ID_B', 'REFSEQ_PROTEIN_B',
             'EVIDENCE', 'PMIDS']


def process_flat_files(id_mappings_file, complexes_file=None, ptm_file=None,
                       seq_file=None, ppi_file=None):
    id_df = pd.read_csv(id_mappings_file, delimiter='\t', names=_hprd_id_cols,
                        dtype='str')
    id_df = id_df.set_index('HPRD_ID')
    if complexes_file is None and  ptm_file is None and ppi_file is None:
        raise ValueError('At least one of complexes_file, ptm_file, or '
                         'ppi_file must be given.')
    if ptm_file and not seq_file:
        raise ValueError('If ptm_file is given, seq_file must also be given.')
    # Load complexes into dataframe
    cplx_df = None
    if complexes_file:
        cplx_df = pd.read_csv(complexes_file, delimiter='\t', names=_cplx_cols,
                              dtype='str')
    # Load ptm data into dataframe
    ptm_df = None
    seq_dict = None
    if ptm_file:
        ptm_df = pd.read_csv(ptm_file, delimiter='\t', names=_ptm_cols,
                             dtype='str', na_values='-')
        # Load protein sequences as a dict keyed by RefSeq ID
        seq_dict = load_fasta_sequences(seq_file, id_index=2)
    # Load the PPI data into dataframe
    ppi_df = None
    if ppi_file:
        ppi_df = pd.read_csv(ppi_file, delimiter='\t', names=_ppi_cols,
                             dtype='str')
    # Create the processor
    return HprdProcessor(id_df, cplx_df, ptm_df, seq_dict, ppi_df)

