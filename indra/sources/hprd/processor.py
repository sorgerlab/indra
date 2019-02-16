import logging
from numpy import nan
from indra.databases import hgnc_client
from indra.statements import *


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

    def __init__(self, id_df, cplx_df=None, ptm_df=None):
        if cplx_df is None and ptm_df is None:
            raise ValueError('At least one of cplx_df or phospho_df must be '
                             'specified.')
        self.statements = []
        self.id_df = id_df
        if cplx_df is not None:
            self.get_complexes(cplx_df)
        if ptm_df is not None:
            self.get_ptms(ptm_df)

    def get_complexes(self, cplx_df):
        # Bring the ID information from the ID table into the dataframe of
        # complexes
        for cplx_id, this_cplx in cplx_df.groupby('CPLX_ID'):
            agents = []
            for hprd_id in this_cplx.HPRD_ID:
                ag = self._make_agent_from_hprd_id(hprd_id)
                if ag is not None:
                    agents.append(ag)
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
            sub_ag = self._make_agent_from_hprd_id(row['HPRD_ID'])
            enz_id = _nan_to_none(row['ENZ_HPRD_ID'])
            enz_ag = self._make_agent_from_hprd_id(enz_id)
            res = _nan_to_none(row['RESIDUE'])
            pos = _nan_to_none(row['POSITION'])
            if pos is not None and ';' in pos:
                pos, dash = pos.split(';')
                assert dash == '-'
            # Sites are semicolon delimited
            ev_list = self._get_evidence(
                    row['HPRD_ID'], row['HPRD_ISOFORM'], row['PMIDS'],
                    row['EVIDENCE'], 'ptms')
            stmt = ptm_class(enz_ag, sub_ag, res, pos, evidence=ev_list)
            self.statements.append(stmt)

    def _make_agent_from_hprd_id(self, hprd_id):
        if hprd_id is None:
            return None
        # Get the Entrez ID from the ID mappings dataframe
        egid = self.id_df.loc[hprd_id].EGID
        refseq_id = self.id_df.loc[hprd_id].REFSEQ_PROTEIN
        # Get the HGNC ID
        hgnc_id = hgnc_client.get_hgnc_from_entrez(egid)
        if not hgnc_id:
            return None
        # Get the (possibly updated) HGNC Symbol
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        assert hgnc_name is not None
        # Get the Uniprot ID, if present (there is no UP id for TRA or TRB)
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        if not up_id:
            logger.info("Skipping entry for Entrez ID %s->HGNC %s "
                        "with no Uniprot ID" % (egid, hgnc_name))
        # Make db_refs, return Agent
        db_refs = {'HGNC': hgnc_id, 'UP': up_id, 'EGID': egid,
                   'REFSEQ_PROT': refseq_id}
        return Agent(hgnc_name, db_refs=db_refs)

    def _make_agent_from_refseq_id(self, refseq_id):
        if refseq_id is None:
            return None
        # Get the Uniprot IDs from the uniprot client
        up_ids = uniprot_client.get_ids_from_refseq(refseq_id,
                                                    reviewed_only=True)
        if len(up_ids) == 0:
            print("No reviewed UP IDs for refseq_id %s" % refseq_id)
        if len(up_ids) > 1:
            print("More than 1 UP ID for refseq_id %s" % refseq_id)
        return None

    def _get_evidence(self, hprd_id, isoform_id, pmid_str, evidence_type,
                      info_type):
        pmids = pmid_str.split(',')
        ev_list = []
        for pmid in pmids:
            ev_type_list = evidence_type.split(';')
            ev = Evidence(source_api='hprd',
                          source_id=_hprd_url(hprd_id, isoform_id, info_type),
                          annotations={'evidence': ev_type_list},
                          pmid=pmid)
            ev_list.append(ev)
        return ev_list


def _hprd_url(hprd_id, isoform_id, info_type):
    if info_type == 'interactions':
        return ('http://hprd.org/interactions?hprd_id=%s&isoform_id=%s_1'
                '&isoform_name=' % (hprd_id, hprd_id))
    else:
        return None


def _nan_to_none(val):
    return None if val is nan else val

