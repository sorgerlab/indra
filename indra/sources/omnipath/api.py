import logging
import requests
from .processor import OmniPathProcessor

logger = logging.getLogger("omnipath")


op_url = 'http://omnipathdb.org'


def process_from_web():
    ptm_json = _get_modifications()
    ligrec_json = _get_interactions()
    return OmniPathProcessor(ptm_json=ptm_json, ligrec_json=ligrec_json)


def _get_modifications():
    """Get all PTMs from Omnipath in JSON format.

    Returns
    -------
    JSON content for PTMs.
    """
    #params = {'format': 'json', 'fields':['sources', 'references']}
    params = {'format': 'json', 'fields':['sources']}
    ptm_url = '%s/ptms' % op_url
    res = requests.get(ptm_url, params=params)
    if not res.status_code == 200 or not res.text:
        return None
    else:
        return res.json()


def _get_interactions(datasets=None):
    """Wrapper for calling the interactions omnipath interactions API

    Some interesting options:
    format: json,tab,table,text,tsv
    datasets: dorothea,kinaseextra,ligrecextra,lncrna_mrna,mirnatarget,
              omnipath,pathwayextra,tf_mirna,tf_target,tfregulons
    types: lncrna_post_transcriptional,mirna_transcriptional,
           post_transcriptional,post_translational,transcriptional
    databases: ABS,ACSN,ARN,Adhesome,AlzPathway,Baccin2019,BioGRID,CA1,
               CancerCellMap,CellPhoneDB,DEPOD,DIP,DOMINO,DeathDomain,
               DoRothEA_A,DoRothEA_B,DoRothEA_C,DoRothEA_D,ELM,EMBRACE,
               ENCODE-distal,ENCODE-proximal,ENCODE_tf-mirna,Guide2Pharma,
               HPMR,HPRD-phos,HTRIdb,ICELLNET,InnateDB,KEA,KEGG,KEGG-MEDICUS,
               Kirouac2010,LMPID,LRdb,Li2012,LncRNADisease,MIMP,MPPI,
               Macrophage,MatrixDB,NRF2ome,NetPath,ORegAnno,PAZAR,PDZBase,
               PhosphoNetworks,PhosphoPoint,PhosphoSite,PhosphoSite_noref,
               ProtMapper,Ramilowski2015,SIGNOR,SPIKE,SignaLink3,TRIP,TransmiR,
               Wang,dbPTM,iPTMnet,iTALK,lncrnadb,miR2Disease,miRDeathDB,
               miRTarBase,miRecords,ncRDeathDB,phosphoELM
    genesymbols: 0,1,no,yes
    fields: curation_effort, databases, dorothea_chipseq, dorothea_coexp,
            dorothea_curated, dorothea_level, dorothea_tfbs, entity_type,
            ncbi_tax_id, organism, references, resources, sources,
            tfregulons_chipseq, tfregulons_coexp, tfregulons_curated,
            tfregulons_level, tfregulons_tfbs, type
    organisms: 10090,10116,9606 (9606 is human)
    directed: 0,1,no,yes
    signed: 0,1,no,yes
    entity_types: complex,lncrna,mirna,protein,small_molecule

    See full list of options here:
    https://omnipathdb.org/queries/interactions

    Parameters
    ----------
    datasets
        A list of dataset names. Options are:
            dorothea, kinaseextra, ligrecextra, lncrna_mrna, mirnatarget,
            omnipath, pathwayextra, tf_mirna, tf_target, tfregulons
        Default: 'ligrecextra'

    Returns
    -------
    dict
        json of database request
    """
    interactions_url = '%s/interactions' % op_url
    params = {
        'fields': ['curation_effort', 'entity_type', 'references',
                   'resources', 'sources', 'type'],
        'format': 'json',
        'datasets': datasets or ['ligrecextra']
    }
    res = requests.get(interactions_url, params=params)
    res.raise_for_status()

    return res.json()
