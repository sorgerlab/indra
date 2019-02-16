import logging
from numpy import nan
from collections import Counter
from protmapper import uniprot_client
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

    def __init__(self, id_df, cplx_df=None, ptm_df=None):
        if cplx_df is None and ptm_df is None:
            raise ValueError('At least one of cplx_df or phospho_df must be '
                             'specified.')
        self.statements = []
        self.id_df = id_df
        # Keep track of the ID mapping issues encountered
        self.no_hgnc_for_egid = []
        self.no_hgnc_for_egid_or_name = []
        self.no_up_for_hgnc = []
        self.no_up_for_refseq = []
        self.many_ups_for_refseq = []

        if cplx_df is not None:
            self.get_complexes(cplx_df)
        if ptm_df is not None:
            self.get_ptms(ptm_df)

        self.no_hgnc_for_egid = Counter(self.no_hgnc_for_egid)
        self.no_hgnc_for_egid_or_name = Counter(self.no_hgnc_for_egid_or_name)
        self.no_up_for_hgnc = Counter(self.no_up_for_hgnc)
        self.no_up_for_refseq = Counter(self.no_up_for_refseq)
        self.many_ups_for_refseq = Counter(self.many_ups_for_refseq)



    def get_complexes(self, cplx_df):
        # Bring the ID information from the ID table into the dataframe of
        # complexes
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
            # Sites are semicolon delimited
            ev_list = self._get_evidence(
                    row['HPRD_ID'], row['HPRD_ISOFORM'], row['PMIDS'],
                    row['EVIDENCE'], 'ptms')
            stmt = ptm_class(enz_ag, sub_ag, res, pos, evidence=ev_list)
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
        # the Entrez ID has been discontinued or replaced. As a hail mary
        # we try to get the HGNC ID from the gene symbol
        if not hgnc_id:
            self.no_hgnc_for_egid.append(egid)
            logger.info('No HGNC ID for HPRD ID %s, EGID %s' % (hprd_id, egid))
            hgnc_symbol = self.id_df.loc[hprd_id].HGNC_SYMBOL
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
            # If we still couldn't get an HGNC ID, then bail out
            if not hgnc_id:
                logger.info('No HGNC ID from EGID %s or from gene symbol %s'
                            % (egid, hgnc_symbol))
                self.no_hgnc_for_egid_or_name.append((egid, hgnc_symbol))
                return None
        # Get the (possibly updated) HGNC Symbol
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        assert hgnc_name is not None
        # If we have provided the RefSeq ID, it's because we need to make
        # sure that we are getting the right isoform-specific ID (for sequence
        # positions of PTMs). In this case we get the Uniprot ID from the
        # Refseq->UP mappings in the protmapper.uniprot_client.
        if refseq_id is not None:
            # Get the Uniprot IDs from the uniprot client
            up_ids = uniprot_client.get_ids_from_refseq(refseq_id,
                                                        reviewed_only=True)
            if len(up_ids) == 0:
                #logger.info("No reviewed UP IDs for refseq_id %s" % refseq_id)
                up_id = None
                self.no_up_for_refseq.append(refseq_id)
            else:
                if len(up_ids) > 1:
                    #logger.info("More than 1 UP ID for refseq_id %s" %
                    #            refseq_id)
                    self.many_ups_for_refseq.append(refseq_id)
                up_id = up_ids[0]

        # On the other hand, if we haven't passed in the Refseq ID, it is
        # because we can obtain it from the HPRD ID table; and we can obtain
        # the Uniprot ID from HGNC.
        else:
            # Get Refseq ID from the HPRD ID table
            refseq_id = self.id_df.loc[hprd_id].REFSEQ_PROTEIN
            # Get Uniprot ID from HGNC, if possible (e.g. there is no UP id
            # for genes TRA or TRB)
            up_id = hgnc_client.get_uniprot_id(hgnc_id)
            if not up_id:
                self.no_up_for_hgnc.append((hgnc_name, hgnc_id))
                logger.info("Skipping entry for Entrez ID %s->HGNC %s "
                            "with no Uniprot ID" % (egid, hgnc_name))
        # Make db_refs, return Agent
        db_refs = {'HGNC': hgnc_id, 'UP': up_id, 'EGID': egid,
                   'REFSEQ_PROT': refseq_id}
        return Agent(hgnc_name, db_refs=db_refs)


    def _get_evidence(self, hprd_id, isoform_id, pmid_str, evidence_type,
                      info_type):
        pmids = pmid_str.split(',')
        ev_list = []
        for pmid in pmids:
            ev_type_list = [] if evidence_type is nan \
                               else evidence_type.split(';')
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

