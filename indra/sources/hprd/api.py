__all__ = ['process_archive', 'process_flat_files']

import requests
import tarfile
import pandas as pd
from protmapper.uniprot_client import load_fasta_sequences, \
    load_fasta_sequence_lines
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


file_mappings = {
    'id_mappings_file': 'HPRD_ID_MAPPINGS.txt',
    'complexes_file': 'PROTEIN_COMPLEXES.txt',
    'ptm_file': 'POST_TRANSLATIONAL_MODIFICATIONS.txt',
    'ppi_file': 'BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt',
    'seq_file': 'PROTEIN_SEQUENCES.txt',
}


def process_archive(fname):
    """Get INDRA Statements from HPRD data in a single tar.gz file.

    The latest release, HPRD_FLAT_FILES_041310.tar.gz can be downloaded from
    http://hprd.org/download after registration.

    Parameters
    ----------
    fname : str
        Path to HPRD tar.gz file.

    Returns
    -------
    HprdProcessor
        An HprdProcessor object which contains a list of extracted INDRA
        Statements in its statements attribute.
    """
    with tarfile.open(fname, "r:gz") as fh:
        prefix = fh.next().name.split('/')[0]
        files = {k: fh.extractfile(prefix + '/' + v) for k, v in
                 file_mappings.items()}
        return process_flat_files(**files)


def process_flat_files(id_mappings_file, complexes_file=None, ptm_file=None,
                       ppi_file=None, seq_file=None, motif_window=7):
    """Get INDRA Statements from HPRD data in individual files.

    Of the arguments, `id_mappings_file` is required, and at least one of
    `complexes_file`, `ptm_file`, and `ppi_file` must also be given.  If
    `ptm_file` is given, `seq_file` must also be given.

    Note that many proteins (> 1,600) in the HPRD content are associated with
    outdated RefSeq IDs that cannot be mapped to Uniprot IDs. For these, the
    Uniprot ID obtained from the HGNC ID (itself obtained from the Entrez ID)
    is used. Because the sequence referenced by the Uniprot ID obtained this
    way may be different from the (outdated) RefSeq sequence included with the
    HPRD content, it is possible that this will lead to invalid site positions
    with respect to the Uniprot IDs.

    To allow these site positions to be mapped during assembly, the
    Modification statements produced by the HprdProcessor include an additional
    key in the `annotations` field of their Evidence object. The annotations
    field is called 'site_motif' and it maps to a dictionary with three
    elements: 'motif', 'respos', and 'off_by_one'. 'motif' gives the peptide
    sequence obtained from the RefSeq sequence included with HPRD. 'respos'
    indicates the position in the peptide sequence containing the residue.
    Note that these positions are ONE-INDEXED (not zero-indexed). Finally, the
    'off-by-one' field contains a boolean value indicating whether the correct
    position was inferred as being an off-by-one (methionine cleavage) error.
    If True, it means that the given residue could not be found in the HPRD
    RefSeq sequence at the given position, but a matching residue was found at
    position+1, suggesting a sequence numbering based on the methionine-cleaved
    sequence. The peptide included in the 'site_motif' dictionary is based on
    this updated position.

    Parameters
    ----------
    id_mappings_file : str
        Path to HPRD_ID_MAPPINGS.txt file.
    complexes_file : Optional[str]
        Path to PROTEIN_COMPLEXES.txt file.
    ptm_file : Optional[str]
        Path to POST_TRANSLATIONAL_MODIFICATIONS.txt file.
    ppi_file : Optional[str]
        Path to BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt file.
    seq_file : Optional[str]
        Path to PROTEIN_SEQUENCES.txt file.
    motif_window : int
        Number of flanking amino acids to include on each side of the
        PTM target residue in the 'site_motif' annotations field of the
        Evidence for Modification Statements. Default is 7.

    Returns
    -------
    HprdProcessor
        An HprdProcessor object which contains a list of extracted INDRA
        Statements in its statements attribute.
    """
    id_df = pd.read_csv(id_mappings_file, delimiter='\t', names=_hprd_id_cols,
                        dtype='str')
    id_df = id_df.set_index('HPRD_ID')
    if complexes_file is None and ptm_file is None and ppi_file is None:
        raise ValueError('At least one of complexes_file, ptm_file, or '
                         'ppi_file must be given.')
    if ptm_file and not seq_file:
        raise ValueError('If ptm_file is given, seq_file must also be given.')
    # Load complexes into dataframe
    cplx_df = None
    if complexes_file:
        cplx_df = pd.read_csv(complexes_file, delimiter='\t', names=_cplx_cols,
                              dtype='str', na_values=['-', 'None'])
    # Load ptm data into dataframe
    ptm_df = None
    seq_dict = None
    if ptm_file:
        ptm_df = pd.read_csv(ptm_file, delimiter='\t', names=_ptm_cols,
                             dtype='str', na_values='-')
        # Load protein sequences as a dict keyed by RefSeq ID
        # If this comes from a tar.gz archive directly, we need to get the
        # lines first and decode them, otherwise we can cal a function
        # that expects a standalone file path
        if hasattr(seq_file, 'readlines'):
            lines = [l.decode('utf-8') for l in seq_file.readlines()]
            seq_dict = load_fasta_sequence_lines(lines, id_index=2)
        else:
            seq_dict = load_fasta_sequences(seq_file, id_index=2)
    # Load the PPI data into dataframe
    ppi_df = None
    if ppi_file:
        ppi_df = pd.read_csv(ppi_file, delimiter='\t', names=_ppi_cols,
                             dtype='str')
    # Create the processor
    return HprdProcessor(id_df, cplx_df, ptm_df, ppi_df, seq_dict, motif_window)
