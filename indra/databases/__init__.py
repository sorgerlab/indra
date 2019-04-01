import logging

logger = logging.getLogger(__name__)


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
    identifiers_url = 'http://identifiers.org/'
    bel_scai_url = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/'
    if db_name == 'UP':
        url = identifiers_url + 'uniprot/%s' % db_id
    elif db_name == 'HGNC':
        url = identifiers_url + 'hgnc/HGNC:%s' % db_id
    elif db_name == 'IP':
        url = identifiers_url + 'interpro/%s' % db_id
    elif db_name == 'IPR':
        url = identifiers_url + 'interpro/%s' % db_id
    elif db_name == 'CHEBI':
        url = identifiers_url + 'chebi/%s' % db_id
    elif db_name == 'NCIT':
        url = identifiers_url + 'ncit/%s' % db_id
    elif db_name == 'GO':
        if db_id.startswith('GO:'):
            url = identifiers_url + 'go/%s' % db_id
        else:
            url = identifiers_url + 'go/GO:%s' % db_id
    elif db_name in ('PUBCHEM', 'PCID'):  # Assuming PCID = PubChem compound ID
        if db_id.startswith('PUBCHEM:'):
            db_id = db_id[8:]
        elif db_id.startswith('PCID:'):
            db_id = db_id[5:]
        url = identifiers_url + 'pubchem.compound/%s' % db_id
    elif db_name == 'PF':
        url = identifiers_url + 'pfam/%s' % db_id
    elif db_name == 'MIRBASEM':
        url = identifiers_url + 'mirbase.mature/%s' % db_id
    elif db_name == 'MIRBASE':
        url = identifiers_url + 'mirbase/%s' % db_id
    elif db_name == 'MESH':
        url = identifiers_url + 'mesh/%s' % db_id
    elif db_name == 'EGID':
        url = identifiers_url + 'ncbigene/%s' % db_id
    elif db_name == 'HMDB':
        url = identifiers_url + 'hmdb/%s' % db_id
    elif db_name == 'LINCS':
        if db_id.startswith('LSM-'):  # Lincs Small Molecule ID
            url = identifiers_url + 'lincs.smallmolecule/%s' % db_id
        elif db_id.startswith('LCL-'):  # Lincs Cell Line ID
            url = identifiers_url + 'lincs.cell/%s' % db_id
        else:  # Assume LINCS Protein
            url = identifiers_url + 'lincs.protein/%s' % db_id
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
    elif db_name == 'FPLX':
        url = 'http://identifiers.org/fplx/%s' % db_id
    elif db_name == 'LNCRNADB':
        if db_id.startswith('ENSG'):
            url = 'http://www.lncrnadb.org/search/?q=%s' % db_id
        else:  # Assmuing HGNC symbol
            url = 'http://www.lncrnadb.org/%s/' % db_id
    elif db_name == 'NXPFA':
        url = 'https://www.nextprot.org/term/FA-%s' % db_id
    elif db_name in ('UN', 'WDI', 'FAO'):
        url = 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/%s' % \
                db_id
    elif db_name == 'HUME':
        url = ('https://github.com/BBN-E/Hume/blob/master/resource/ontologies/'
               'hume_ontology/%s' % db_id)
    elif db_name == 'CWMS':
        url = 'http://trips.ihmc.us/%s' % db_id
    elif db_name == 'SIGNOR':  # Assuming db_id == Primary ID
        url = 'https://signor.uniroma2.it/relation_result.php?id=%s' % db_id
    elif db_name == 'SOFIA':
        url = 'http://cs.cmu.edu/sofia/%s' % db_id
    elif db_name == 'CHEMBL':
        if not db_id.startswith('CHEMBL'):
            db_id = 'CHEMBL%s' % db_id
        url = identifiers_url + 'chembl.compound/%s' % db_id
    elif db_name == 'NONCODE':
        url = 'http://www.noncode.org/show_gene.php?id=NONHSAG%s' % db_id
    elif db_name == 'TEXT':
        return None
    else:
        logger.warning('Unhandled name space %s' % db_name)
        url = None
    return url
