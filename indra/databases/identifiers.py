import re
import logging
from indra.resources import load_resource_json


logger = logging.getLogger(__name__)


identifiers_url = 'https://identifiers.org'

# These are just special cases of name spaces where the mapping from INDRA to
# identifiers.org is not a question of simplecapitalization.
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
    'CVCL': 'cellosaurus',
}

# These are namespaces used by INDRA that don't have corresponding
# identifiers.org entries
non_registry = {
    'SDIS', 'SCHEM', 'SFAM', 'SCOMP', 'HMS-LINCS', 'NXPFA',
    'OMIM', 'LSPCI', 'UPLOC', 'BFO', 'CCLE', 'CLO', 'GENBANK',
    'DRUGBANK.SALT', 'SMILES',
}

# These are reverse mappings from identifiers.org namespaces to INDRA
# namespaces
identifiers_reverse = {
    v: k for k, v in identifiers_mappings.items()
}
# We have to patch this one because it is ambiguous
identifiers_reverse['ncbigene'] = 'EGID'

# These are only the URLs that are strictly prefixes and not more complicated
# patterns. This is because some downstream code uses these as prefixes
# rather than arbitrary patterns.
url_prefixes = {
    # Biology namespaces
    'NXPFA': 'https://www.nextprot.org/term/FA-',
    'SIGNOR': 'https://signor.uniroma2.it/relation_result.php?id=',
    'LSPCI': 'https://labsyspharm.github.io/lspci/',
    # WM namespaces
    'UN': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'WDI': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'FAO': 'https://github.com/clulab/eidos/wiki/JSON-LD#Grounding/',
    'HUME': ('https://github.com/BBN-E/Hume/blob/master/resource/ontologies'
             '/hume_ontology/'),
    'CWMS': 'http://trips.ihmc.us/',
    'SOFIA': 'http://cs.cmu.edu/sofia/',
}


def get_ns_from_identifiers(identifiers_ns):
    """"Return a namespace compatible with INDRA from an identifiers namespace.

    For example, this function can be used to map 'uniprot' to 'UP'.

    Parameters
    ----------
    identifiers_ns : str
        An identifiers.org standard namespace.

    Returns
    -------
    str or None
        The namespace compatible with INDRA's internal representation or
        None if the given namespace isn't an identifiers.org standard.
    """
    reg_entry = identifiers_registry.get(identifiers_ns.lower())
    if not reg_entry:
        return None
    mapping = identifiers_reverse.get(identifiers_ns.lower())
    if mapping:
        return mapping
    else:
        return identifiers_ns.upper()


def get_ns_id_from_identifiers(identifiers_ns, identifiers_id):
    """Return a namespace/ID pair compatible with INDRA from identifiers.

    Parameters
    ----------
    identifiers_ns : str
        An identifiers.org standard namespace.
    identifiers_id : str
        An identifiers.org standard ID in the given namespace.

    Returns
    -------
    (str, str)
        A namespace and ID that are valid in INDRA db_refs.
    """
    reg_entry = identifiers_registry.get(identifiers_ns.lower())
    db_ns = get_ns_from_identifiers(identifiers_ns)
    if db_ns is None:
        return None, None
    db_id = identifiers_id
    if reg_entry['namespace_embedded']:
        if not identifiers_id.startswith(identifiers_ns.upper()):
            db_id = '%s:%s' % (identifiers_ns.upper(), identifiers_id)
    return db_ns, db_id


def get_identifiers_ns(db_name):
    """Map an INDRA namespace to an identifiers.org namespace when possible.

    Example: this can be used to map 'UP' to 'uniprot'.

    Parameters
    ----------
    db_name : str
        An INDRA namespace to map to identifiers.org

    Returns
    -------
    str or None
        An identifiers.org namespace or None if not available.
    """
    mapped_db_name = identifiers_mappings.get(db_name, db_name.lower())
    if mapped_db_name not in identifiers_registry:
        return None
    return mapped_db_name


def get_url_prefix(db_name):
    """Return the URL prefix for a given namespace."""
    identifiers_ns = get_identifiers_ns(db_name)

    if identifiers_ns:
        identifiers_entry = identifiers_registry.get(identifiers_ns)
        if not identifiers_entry['namespace_embedded']:
            return '%s/%s:' % (identifiers_url, identifiers_ns.lower())
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


def namespace_embedded(db_ns: str) -> bool:
    """Return true if this namespace requires IDs to have namespace embedded.

    This function first maps the given namespace to an identifiers.org
    namespace and then checks the registry to see if namespaces need
    to be embedded in IDs. If yes, it returns True. If not, or the ID can't
    be mapped to identifiers.org, it returns False

    Parameters
    ----------
    db_ns :
        The namespace to check.

    Returns
    -------
    :
        True if the namespace is known to be embedded in IDs of this
        namespace. False if unknown or known not to be embedded.
    """
    identifiers_ns = get_identifiers_ns(db_ns)
    if identifiers_ns:
        identifiers_entry = identifiers_registry.get(identifiers_ns)
        if identifiers_entry['namespace_embedded']:
            return True
    return False


def ensure_prefix_if_needed(db_ns: str, db_id: str) -> str:
    """Return an ID ensuring a namespace prefix if known to be needed.

    Parameters
    ----------
    db_ns :
        The namespace associated with the identifier.
    db_id :
        The original identifier.

    Returns
    -------
    :
        The identifier with namespace embedded if needed.
    """
    if namespace_embedded(db_ns):
        return ensure_prefix(db_ns, db_id)
    return db_id


def ensure_prefix(db_ns, db_id, with_colon=True):
    """Return a valid ID that has the given namespace embedded.

    This is useful for namespaces such as CHEBI, GO or BTO that require
    the namespace to be part of the ID. Note that this function always
    ensures that the given db_ns is embedded in the ID and can handle the
    case whene it's already present.

    Parameters
    ----------
    db_ns : str
        A namespace.
    db_id : str
        An ID within that namespace which should have the namespace
        as a prefix in it.
    with_colon: Optional[bool]
        If True, the namespace prefix is followed by a colon in the ID (e.g.,
        CHEBI:12345). Otherwise, no colon is added (e.g., CHEMBL1234).
        Default: True
    """
    if db_id is None:
        return None
    colon = ':' if with_colon else ''
    if not db_id.startswith(f'{db_ns}{colon}'):
        return f'{db_ns}{colon}{db_id}'
    return db_id


def ensure_chebi_prefix(chebi_id):
    """Return a valid CHEBI ID that has the appropriate CHEBI: prefix."""
    return ensure_prefix('CHEBI', chebi_id)


def ensure_chembl_prefix(chembl_id):
    """Return a valid CHEMBL ID that has the appropriate CHEMBL prefix."""
    return ensure_prefix('CHEMBL', chembl_id, with_colon=False)


def _load_identifiers_registry():
    identifiers_registry = load_resource_json('identifiers_patterns.json')
    # Override pattern otherwise patterns like 1.1 can't be used
    # TODO should we allow identifiers like 6.3.2.n3?
    identifiers_registry['ec-code']['pattern'] = \
        '^\\d{1,2}(\\.[n]?\\d{0,3}){0,3}$'
    identifiers_registry['mondo'] = {
        "pattern": "^\\d+$",
        "namespace_embedded": False,
    }
    for value in identifiers_registry.values():
        value["pattern_compiled"] = re.compile(value["pattern"])
    return identifiers_registry


identifiers_registry = _load_identifiers_registry()
