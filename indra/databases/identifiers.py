import re
import json
import logging
from urllib import parse
from indra.resources import load_resource_json


logger = logging.getLogger(__name__)


identifiers_url = 'https://identifiers.org'

identifiers_mappings = {
    'UP': 'uniprot',
    'UPPRO': 'uniprot.chain',
    'UPISO': 'uniprot.isoform',
    'REFSEQ_PROT': 'refseq',
    'PF': 'pfam',
    'IP': 'interpro',
    'ECCODE': 'ec-code',
    'NONCODE': 'noncodev4.rna',
    'LNCRNADB': 'rnacentral',
    'MIRBASEM': 'mirbase.mature',
    'EGID': 'ncbigene',
    'NCBI': 'ncibgene',
    'HGNC_GROUP': 'hgnc.genefamily',
    'LINCS': 'lincs.smallmolecule',
    'PUBCHEM': 'pubchem.compound',
    'CHEMBL': 'chembl.compound',
    'CTD': 'ctd.chemical',
}

# Get reverse mappings and patch one entry to make it unique
identifiers_reverse = {
    v: k for k, v in identifiers_mappings.items()
}

identifiers_reverse['ncbigene'] = 'EGID'

non_registry = {
    'SDIS', 'SCHEM', 'SFAM', 'SCOMP', 'SIGNOR', 'HMS-LINCS', 'NXPFA',
    'OMIM', 'LSPCI', 'UPLOC'
}

non_grounding = {
    'TEXT', 'TEXT_NORM'
}


# These are only the URLs that are strictly prefixes and not more complicated
# patterns. This is because some downstream code uses these as prefixes
# rather than arbitrary patterns.
url_prefixes = {
    # Biology namespaces
    'NXPFA': 'https://www.nextprot.org/term/FA-',
    'SIGNOR': 'https://signor.uniroma2.it/relation_result.php?id=',
    # WM namespaces
    'UN': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'WDI': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'FAO': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'HUME': ('https://github.com/BBN-E/Hume/blob/master/resource/ontologies'
             '/hume_ontology/'),
    'CWMS': 'http://trips.ihmc.us/',
    'SOFIA': 'http://cs.cmu.edu/sofia/',
}


def get_ns_id_from_identifiers(identifiers_ns, identifiers_id):
    reg_entry = identifiers_registry.get(identifiers_ns.lower())
    if not reg_entry:
        return None, None
    mapping = identifiers_reverse.get(identifiers_ns.lower())
    if mapping:
        db_ns = mapping
    else:
        db_ns = identifiers_ns.upper()
    db_id = identifiers_id
    if reg_entry['namespace_embedded']:
        if not identifiers_id.startswith(identifiers_ns.upper()):
            db_id = '%s:%s' % (identifiers_ns.upper(), identifiers_id)
    return db_ns, db_id


def get_url_prefix(db_name):
    mapped_db_name = identifiers_mappings.get(db_name, db_name.lower())
    identifiers_entry = identifiers_registry.get(mapped_db_name)
    if identifiers_entry:
        if not identifiers_entry['namespace_embedded']:
            return '%s/%s:' % (identifiers_url, mapped_db_name.lower())
        else:
            return '%s/' % identifiers_url
    else:
        if db_name in url_prefixes:
            return url_prefixes[db_name]
    return None


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
    if db_name == 'CHEMBL':
        db_id = ensure_chembl_prefix(db_id)
    elif db_name == 'CHEBI':
        db_id = ensure_chebi_prefix(db_id)

    prefix = get_url_prefix(db_name)
    if prefix:
        return '%s%s' % (prefix, db_id)

    # Otherwise, we have to handle some special cases
    bel_scai_url = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/'
    if db_name == 'LINCS':
        if db_id.startswith('LSM-'):  # Lincs Small Molecule ID
            url = identifiers_url + '/lincs.smallmolecule:%s' % db_id
        elif db_id.startswith('LCL-'):  # Lincs Cell Line ID
            url = identifiers_url + '/lincs.cell:%s' % db_id
        else:  # Assume LINCS Protein
            url = identifiers_url + '/lincs.protein:%s' % db_id
    elif db_name == 'CHEMBL':
        if not db_id.startswith('CHEMBL'):
            db_id = 'CHEMBL%s' % db_id
        url = identifiers_url + '/chembl.compound:%s' % db_id
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
    elif db_name == 'TEXT' or db_name == 'TEXT_NORM':
        return None
    else:
        logger.warning('Unhandled name space %s' % db_name)
        url = None
    return url


def parse_identifiers_url(url):
    """Retrieve database name and ID given the URL.

    Parameters
    ----------
    url : str
        An identifiers.org URL to parse.

    Returns
    -------
    db_name : str
        An internal database name: HGNC, UP, CHEBI, etc. corresponding to the
        given URL.
    db_id : str
        An identifier in the database.
    """
    # Try matching by string pattern
    db_ns, db_id = None, None
    url_pattern = \
        r'(?:https?)://identifiers.org/([A-Za-z0-9.-]+)(/|:)([A-Za-z0-9:_.-]+)'
    match = re.match(url_pattern, url)
    if match is not None:
        g = match.groups()
        if len(g) == 3:
            pattern_ns, pattern_id = g[0], g[2]
            db_ns, db_id = get_ns_id_from_identifiers(pattern_ns, pattern_id)
            if db_ns == 'HGNC':
                if db_id.startswith('HGNC:'):
                    db_id = db_id[5:]
            # If we got UP and UPPRO, return UPPRO
            if db_ns == 'UP' and '#PRO_' in url:
                db_ns = 'UPPRO'
                db_id = url[url.find('PRO_'):]
            if db_ns and db_id:
                return db_ns, db_id
    for ns, prefix in url_prefixes.items():
        if url.startswith(prefix):
            return ns, url[len(prefix):]

    # Handle other special cases
    for part in ['/lincs.smallmolecule', '/lincs.cell', '/lincs.protein']:
        if part in url:
            return 'LINCS', url[(url.find(part) + len(part) + 1):]
    if '/chembl.compound' in url:
        return 'CHEMBL', url[
            (url.find('/chembl.compound') + len('/chembl.compound:')):]
    if 'lincs.hms.harvard.edu' in url:
        return 'HMS-LINCS', url[len('http://lincs.hms.harvard.edu/db/sm/'):-4]
    if 'selventa-legacy-chemicals/' in url:
        return 'SCHEM', None
    if 'selventa-named-complexes/' in url:
        return 'SCOMP', None
    if 'selventa-protein-families/' in url:
        return 'SFAM', None
    else:
        logger.warning('Could not parse URL %s' % url)
    return None, None


def ensure_prefix(db_ns, db_id, with_colon=True):
    if db_id is None:
        return None
    colon = ':' if with_colon else ''
    if not db_id.startswith(f'{db_ns}{colon}'):
        return f'{db_ns}{colon}{db_id}'
    return db_id


def ensure_chebi_prefix(chebi_id):
    return ensure_prefix('CHEBI', chebi_id)


def ensure_chembl_prefix(chembl_id):
    return ensure_prefix('CHEMBL', chembl_id, with_colon=False)


identifiers_registry = load_resource_json('identifiers_patterns.json')
