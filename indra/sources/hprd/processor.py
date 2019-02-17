import logging
from numpy import nan
from collections import Counter
from protmapper import uniprot_client, ProtMapper
from indra.statements import *
from indra.databases import hgnc_client


logger = logging.getLogger('hprd')


# Map of HPRD PTM types to INDRA Statement types
_ptm_map = {
    'ADP-Ribosylation': None,
    'Acetylation': Acetylation,
    'Alkylation': None,
    'Amidation': None,
    'Carboxylation': None,
    'Deacetylation': Deacetylation,
    'Deacylation': None,
    'Deglycosylation': Deglycosylation,
    'Demethylation': Demethylation,
    'Dephosphorylation': Dephosphorylation,
    'Disulfide Bridge': None,
    'Farnesylation': Farnesylation,
    'Glucuronosyl transfer': None,
    'Glutathionylation': None,
    'Glycation': None,
    'Glycosyl phosphatidyl inositol GPI': None,
    'Glycosylation': Glycosylation,
    'Hydroxylation': Hydroxylation,
    'Methylation': Methylation,
    'Myristoylation': Myristoylation,
    'Neddylation': None,
    'Nitration': None,
    'Palmitoylation': Palmitoylation,
    'Phosphorylation': Phosphorylation,
    'Prenylation': None,
    'Proteolytic Cleavage': None,
    'S-Nitrosylation': None,
    'Sulfation': None,
    'Sumoylation': Sumoylation,
    'Transglutamination': None,
    'Ubiquitination': Ubiquitination,
}


class HprdProcessor(object):
    """Get INDRA Statements from HPRD data row objects.

    Parameters
    ----------

    Attributes
    ----------
    statements : list of INDRA Statements
    """

    def __init__(self, id_df, cplx_df=None, ptm_df=None, seq_dict=None,
                 ppi_df=None):
        if cplx_df is None and ptm_df is None and ppi_df is None:
            raise ValueError('At least one of cplx_df, ptm_df, or ppi_df must '
                             'be specified.')
        if ptm_df is not None and not seq_dict:
            raise ValueError('If ptm_df is given, seq_dict must also be given.')

        self.statements = []
        self.id_df = id_df
        self.seq_dict = seq_dict

        # Keep track of the ID mapping issues encountered
        self.no_hgnc_for_egid = []
        self.no_up_for_hgnc = []
        self.no_up_for_refseq = []
        self.many_ups_for_refseq = []
        self.invalid_site_pos = []
        self.off_by_one = []

        # Do the actual processing
        if cplx_df is not None:
            self.get_complexes(cplx_df)
        if ptm_df is not None:
            self.get_ptms(ptm_df)
        if ppi_df is not None:
            self.get_ppis(ppi_df)

        # Tabulate IDs causing issues
        self.no_hgnc_for_egid = Counter(self.no_hgnc_for_egid)
        self.no_up_for_hgnc = Counter(self.no_up_for_hgnc)
        self.no_up_for_refseq = Counter(self.no_up_for_refseq)
        self.many_ups_for_refseq = Counter(self.many_ups_for_refseq)

    def get_complexes(self, cplx_df):
        # Group the agents for the complex
        for cplx_id, this_cplx in cplx_df.groupby('CPLX_ID'):
            agents = []
            for hprd_id in this_cplx.HPRD_ID:
                ag = self._make_agent(hprd_id)
                if ag is not None:
                    agents.append(ag)
            # Make sure we got some agents!
            if not agents:
                continue
            # Get evidence info from first member of complex
            row0 = this_cplx.iloc[0]
            isoform_id = '%s_1' % row0.HPRD_ID
            ev_list = self._get_evidence(row0.HPRD_ID, isoform_id, row0.PMIDS,
                                         row0.EVIDENCE, 'interactions')
            stmt = Complex(agents, evidence=ev_list)
            self.statements.append(stmt)

    def get_ptms(self, ptm_df):
        # Iterate over the rows of the dataframe
        for ix, row in ptm_df.iterrows():
            # Check the modification type; if we can't make an INDRA statement
            # for it, then skip it
            ptm_class = _ptm_map[row['MOD_TYPE']]
            if ptm_class is None:
                continue
            # Use the Refseq protein ID for the substrate to make sure that
            # we get the right Uniprot ID for the isoform
            sub_ag = self._make_agent(row['HPRD_ID'],
                                      refseq_id=row['REFSEQ_PROTEIN'])

            # If we couldn't get the substrate, skip the statement
            if sub_ag is None:
                continue
            enz_id = _nan_to_none(row['ENZ_HPRD_ID'])
            enz_ag = self._make_agent(enz_id)
            res = _nan_to_none(row['RESIDUE'])
            pos = _nan_to_none(row['POSITION'])
            if pos is not None and ';' in pos:
                pos, dash = pos.split(';')
                assert dash == '-'
            # As a fallback for later site mapping, we also get the protein
            # sequence information in case there was a problem with the
            # RefSeq->Uniprot mapping
            assert res
            assert pos
            motif_dict = self._get_seq_motif(row['REFSEQ_PROTEIN'], res, pos)
            # Get evidence
            ev_list = self._get_evidence(
                    row['HPRD_ID'], row['HPRD_ISOFORM'], row['PMIDS'],
                    row['EVIDENCE'], 'ptms', motif_dict)
            stmt = ptm_class(enz_ag, sub_ag, res, pos, evidence=ev_list)
            self.statements.append(stmt)

    def get_ppis(self, ppi_df):
        for ix, row in ppi_df.iterrows():
            agA = self._make_agent(row['HPRD_ID_A'])
            agB = self._make_agent(row['HPRD_ID_B'])
            # If don't get valid agents for both, skip this PPI
            if agA is None or agB is None:
                continue
            isoform_id = '%s_1' % row['HPRD_ID_A']
            ev_list = self._get_evidence(
                    row['HPRD_ID_A'], isoform_id, row['PMIDS'],
                    row['EVIDENCE'], 'interactions')
            stmt = Complex([agA, agB], evidence=ev_list)
            self.statements.append(stmt)

    def _make_agent(self, hprd_id, refseq_id=None):
        if hprd_id is None:
            return None
        # Get the basic info (HGNC name/symbol, Entrez ID) from the
        # ID mappings dataframe
        egid = self.id_df.loc[hprd_id].EGID
        if not egid:
            logger.info('No Entrez ID for HPRD ID %s' % hprd_id)
            return None
        # Get the HGNC ID
        hgnc_id = hgnc_client.get_hgnc_from_entrez(egid)
        # If we couldn't get an HGNC ID for the Entrez ID, this means that
        # the Entrez ID has been discontinued or replaced.
        if not hgnc_id:
            self.no_hgnc_for_egid.append(egid)
            return None
        # Get the (possibly updated) HGNC Symbol
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        assert hgnc_name is not None
        # See if we can get a Uniprot ID from the HGNC symbol--if there is
        # a RefSeq ID we wil also try to use it to get an isoform specific
        # UP ID, but we will have this one to fall back on. But if we can't
        # get one here, then we skip the Statement
        up_id_from_hgnc = hgnc_client.get_uniprot_id(hgnc_id)
        if not up_id_from_hgnc:
            self.no_up_for_hgnc.append((egid, hgnc_name, hgnc_id))
            return None
        # If we have provided the RefSeq ID, it's because we need to make
        # sure that we are getting the right isoform-specific ID (for sequence
        # positions of PTMs). Here we try to get the Uniprot ID from the
        # Refseq->UP mappings in the protmapper.uniprot_client.
        if refseq_id is not None:
            # Get the Uniprot IDs from the uniprot client
            up_ids = uniprot_client.get_ids_from_refseq(refseq_id,
                                                        reviewed_only=True)
            # Nothing for this RefSeq ID (quite likely because the RefSeq ID
            # is obsolete; take the UP ID from HGNC
            if len(up_ids) == 0:
                self.no_up_for_refseq.append(refseq_id)
                up_id = up_id_from_hgnc
            # More than one reviewed entry--no thanks, we'll take the one from
            # HGNC instead
            elif len(up_ids) > 1:
                self.many_ups_for_refseq.append(refseq_id)
                up_id = up_id_from_hgnc
            # We got a unique, reviewed UP entry for the RefSeq ID
            else:
                up_id = up_ids[0]
                # If it's the canonical isoform, strip off the '-1'
                if up_id.endswith('-1'):
                    up_id = up_id.split('-')[0]
        # For completeness, get the Refseq ID from the HPRD ID table
        else:
            refseq_id = self.id_df.loc[hprd_id].REFSEQ_PROTEIN
            up_id = up_id_from_hgnc
        # Make db_refs, return Agent
        db_refs = {'HGNC': hgnc_id, 'UP': up_id, 'EGID': egid,
                   'REFSEQ_PROT': refseq_id}
        return Agent(hgnc_name, db_refs=db_refs)


    def _get_evidence(self, hprd_id, isoform_id, pmid_str, evidence_type,
                      info_type, motif_dict=None):
        pmids = pmid_str.split(',')
        ev_list = []
        if motif_dict is None:
            motif_dict = {}
        for pmid in pmids:
            ev_type_list = [] if evidence_type is nan \
                               else evidence_type.split(';')
            # Add the annotations with the site motif info if relevant
            annotations = {'evidence': ev_type_list}
            annotations.update(motif_dict)
            ev = Evidence(source_api='hprd',
                          source_id=_hprd_url(hprd_id, isoform_id, info_type),
                          pmid=pmid, annotations=annotations)
            ev_list.append(ev)
        return ev_list


    def _get_seq_motif(self, refseq_id, residue, position, window=7):
        seq = self.seq_dict[refseq_id]
        pos_ix = int(position) - 1
        if seq[pos_ix] != residue:
            self.invalid_site_pos.append((refseq_id, residue, position))
            if seq[pos_ix + 1] == residue:
                self.off_by_one.append((refseq_id, residue, position))
            return {}
        else:
            # The index of the residue at the start of the window
            start_ix = pos_ix - window
            motif, respos = ProtMapper.motif_from_position_seq(seq, position)
            return {'site_motif': {'motif': motif, 'respos': respos}}


def _hprd_url(hprd_id, isoform_id, info_type):
    if info_type == 'interactions':
        return ('http://hprd.org/interactions?hprd_id=%s&isoform_id=%s_1'
                '&isoform_name=' % (hprd_id, hprd_id))
    else:
        return None


def _nan_to_none(val):
    return None if val is nan else val

