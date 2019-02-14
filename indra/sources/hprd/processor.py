import logging
from indra.databases import hgnc_client
from indra.statements import Complex, Agent, Evidence


logger = logging.getLogger('hprd')


class HprdProcessor(object):
    """Get INDRA Statements from HPRD data row objects.

    Parameters
    ----------

    Attributes
    ----------
    statements : list of INDRA Statements
    """

    def __init__(self, id_df, cplx_df=None, phospho_df=None):
        if cplx_df is None and phospho_df is None:
            raise ValueError('At least one of cplx_df or phospho_df must be '
                             'specified.')
        self.statements = []
        self.id_df = id_df
        if cplx_df is not None:
            self.get_complexes(cplx_df)


    def get_complexes(self, cplx_df):
        # Bring the ID information from the ID table into the dataframe of
        # complexes
        cplx_id_join = cplx_df.join(self.id_df, 'HPRD_ID', how='left',
                                    lsuffix='CPLX') 
        for cplx_id, this_cplx in cplx_id_join.groupby('CPLX_ID'):
            agents = []
            evidence = []
            for egid in this_cplx.EGID:
                ag = _make_agent_from_egid(egid)
                if ag is not None:
                    agents.append(ag)
            # Get info from first member of complex
            row0 = this_cplx.iloc[0]
            pmids = row0.PMIDS.split(';')
            ev_types = row0.EVIDENCE
            ev_list = []
            for pmid in pmids:
                ev = Evidence(source_api='hprd',
                              source_id=_hprd_url(row0.HPRD_ID, 'interactions'),
                              pmid=pmid)
                ev_list.append(ev)
            stmt = Complex(agents, evidence=ev_list)
            self.statements.append(stmt)


def _hprd_url(hprd_id, info_type):
    if info_type == 'interactions':
        return ('http://hprd.org/interactions?hprd_id=%s&isoform_id=%s_1'
                '&isoform_name=' % (hprd_id, hprd_id))


def _make_agent_from_egid(egid):
    # Make sure the EGID is a string
    egid = str(egid)
    hgnc_id = hgnc_client.get_hgnc_from_entrez(egid)
    if not hgnc_id:
        return None
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert hgnc_name is not None
    up_id = hgnc_client.get_uniprot_id(hgnc_id)
    if not up_id:
        logger.info("Skipping entry for gene %s with no Uniprot ID" %
                    hgnc_name)
    db_refs = {'HGNC': hgnc_id, 'UP': up_id, 'EGID': egid}
    return Agent(hgnc_name, db_refs=db_refs)

