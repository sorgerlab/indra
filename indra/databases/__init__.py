import logging

logger = logging.getLogger('databases')

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
    if db_name == 'UP':
        url = identifiers_url + 'uniprot/%s' % db_id
    elif db_name == 'HGNC':
        url = identifiers_url + 'hgnc/HGNC:%s' % db_id
    elif db_name == 'IP':
        url = identifiers_url + 'interpro/%s' % db_id
    elif db_name == 'CHEBI':
        url = identifiers_url + 'chebi/%s' % db_id
    elif db_name == 'NCIT':
        url = identifiers_url + 'ncit/%s' % db_id
    elif db_name == 'GO':
        url = identifiers_url + 'go/%s' % db_id
    elif db_name == 'PUBCHEM':
        if db_id.startswith('PUBCHEM:'):
            db_id = db_id[8:]
        url = identifiers_url + 'pubchem.compound/%s' % db_id
    elif db_name == 'PF':
        url = identifiers_url + 'pfam/%s' % db_id
    elif db_name == 'MESH':
        url = identifiers_url + 'mesh/%s' % db_id
    elif db_name == 'HMDB':
        url = identifiers_url + 'hmdb/%s' % db_id
    # Special cases with no identifiers entry
    elif db_name == 'FPLX':
        url = 'http://sorger.med.harvard.edu/indra/entities/%s' % db_id
    elif db_name == 'NXPFA':
        url = 'https://www.nextprot.org/term/FA-%s' % db_id
    elif db_name == 'TEXT':
        return None
    else:
        logger.warning('Unhandled name space %s' % db_name)
        url = None
    return url
