import logging

logger = logging.getLogger(__name__)


identifiers_url = 'https://identifiers.org'

# These are only the URLs that are strictly prefixes and not more complicated
# patterns. This is because some downstream code uses these as prefixes
# rather than arbitrary patterns.
url_prefixes = {
    'UP': '%s/uniprot/' % identifiers_url,
    'HGNC': '%s/hgnc/HGNC:' % identifiers_url,
    'IP': '%s/interpro/' % identifiers_url,
    'IPR': '%s/interpro/' % identifiers_url,
    'CHEBI': '%s/chebi/' % identifiers_url,
    'NCIT': '%s/ncit/' % identifiers_url,
    'GO': '%s/go/' % identifiers_url,
    'PCID': '%s/pubchem.compound/' % identifiers_url,
    'PUBCHEM': '%s/pubchem.compound/' % identifiers_url,
    'PF': '%s/pfam/' % identifiers_url,
    'MIRBASEM': '%s/mirbase.mature/' % identifiers_url,
    'MIRBASE': '%s/mirbase/' % identifiers_url,
    'MESH': '%s/mesh/' % identifiers_url,
    'EGID': '%s/ncbigene/' % identifiers_url,
    'HMDB': '%s/hmdb/' % identifiers_url,
    'FPLX': '%s/fplx/' % identifiers_url,
    'REFSEQ_PROT': '%s/refseq:' % identifiers_url,
    'EFO': '%s/efo/' % identifiers_url,
    'HP': '%s/hp/' % identifiers_url,
    'DOID': '%s/' % identifiers_url,  # note that IDs start with DOID:
    'NXPFA': 'https://www.nextprot.org/term/FA-',
    'SIGNOR': 'https://signor.uniroma2.it/relation_result.php?id=',
    'NONCODE': 'http://www.noncode.org/show_gene.php?id=NONHSAG',
    'UN': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'WDI': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'FAO': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'HUME': ('https://github.com/BBN-E/Hume/blob/master/resource/ontologies'
             '/hume_ontology/'),
    'CWMS': 'http://trips.ihmc.us/',
    'SOFIA': 'http://cs.cmu.edu/sofia/',
}


def get_identifiers_url(db_name, db_id):
    """Return an identifiers.org URL for a given database name and ID.

    Parameters
    ----------
    db_name : str
        An internal database name: HGNC, UP, CHEBI, etc.
    db_id : str
        An identifier in the given database.

    Returns
    -------
    url : str
        An identifiers.org URL corresponding to the given database name and ID.
    """
    # This is the case where we have a prefix that we can simply attach the
    # db_id to to get the desired URL.
    prefix = url_prefixes.get(db_name)
    if prefix:
        return '%s%s' % (prefix, db_id)

    # Otherwise, we have to handle some special cases
    bel_scai_url = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/'
    if db_name == 'LINCS':
        if db_id.startswith('LSM-'):  # Lincs Small Molecule ID
            url = identifiers_url + '/lincs.smallmolecule/%s' % db_id
        elif db_id.startswith('LCL-'):  # Lincs Cell Line ID
            url = identifiers_url + '/lincs.cell/%s' % db_id
        else:  # Assume LINCS Protein
            url = identifiers_url + '/lincs.protein/%s' % db_id
    elif db_name == 'CHEMBL':
        if not db_id.startswith('CHEMBL'):
            db_id = 'CHEMBL%s' % db_id
        url = identifiers_url + '/chembl.compound/%s' % db_id
    elif db_name == 'HMS-LINCS':
        url = 'http://lincs.hms.harvard.edu/db/sm/%s-101' % db_id
    # Special cases with no identifiers entry
    elif db_name == 'SCHEM':
        url = bel_scai_url + 'selventa-legacy-chemicals/' + \
            'selventa-legacy-chemicals-20150601.belns'
    elif db_name == 'SCOMP':
        url = bel_scai_url + 'selventa-named-complexes/' + \
            'selventa-named-complexes-20150601.belns'
    elif db_name == 'SFAM':
        url = bel_scai_url + 'selventa-protein-families/' + \
            'selventa-protein-families-20150601.belns'
    elif db_name == 'LNCRNADB':
        if db_id.startswith('ENSG'):
            url = 'http://www.lncrnadb.org/search/?q=%s' % db_id
        else:  # Assmuing HGNC symbol
            url = 'http://www.lncrnadb.org/%s/' % db_id
    elif db_name == 'TEXT':
        return None
    else:
        logger.warning('Unhandled name space %s' % db_name)
        url = None
    return url
