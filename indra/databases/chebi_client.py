import os
import logging
import requests
from lxml import etree
from functools import lru_cache, cmp_to_key
from indra.util import read_unicode_csv
from indra.databases.obo_client import OboClient

_obo_client = OboClient(prefix='chebi')

logger = logging.getLogger(__name__)

# Namespaces used in the XML
chebi_xml_ns = {'n': 'http://schemas.xmlsoap.org/soap/envelope/',
                'c': 'https://www.ebi.ac.uk/webservices/chebi'}


def _add_prefix(chid):
    if chid and not chid.startswith('CHEBI:'):
        return 'CHEBI:%s' % chid
    else:
        return chid


def get_pubchem_id(chebi_id):
    """Return the PubChem ID corresponding to a given ChEBI ID.

    Parameters
    ----------
    chebi_id : str
        ChEBI ID to be converted.

    Returns
    -------
    pubchem_id : str
        PubChem ID corresponding to the given ChEBI ID. If the lookup fails,
        None is returned.
    """
    pubchem_id = chebi_pubchem.get(_add_prefix(chebi_id))
    return pubchem_id


def get_chebi_id_from_pubchem(pubchem_id):
    """Return the ChEBI ID corresponding to a given Pubchem ID.

    Parameters
    ----------
    pubchem_id : str
        Pubchem ID to be converted.

    Returns
    -------
    chebi_id : str
        ChEBI ID corresponding to the given Pubchem ID. If the lookup fails,
        None is returned.
    """
    chebi_id = pubchem_chebi.get(pubchem_id)
    return chebi_id


def get_chembl_id(chebi_id):
    """Return a ChEMBL ID from a given ChEBI ID.

    Parameters
    ----------
    chebi_id : str
        ChEBI ID to be converted.

    Returns
    -------
    chembl_id : str
        ChEMBL ID corresponding to the given ChEBI ID. If the lookup fails,
        None is returned.
    """
    return chebi_chembl.get(_add_prefix(chebi_id))


def get_chebi_id_from_chembl(chembl_id):
    """Return a ChEBI ID from a given ChEBML ID.

    Parameters
    ----------
    chembl_id : str
        ChEBML ID to be converted.

    Returns
    -------
    chebi_id : str
        ChEBI ID corresponding to the given ChEBML ID. If the lookup fails,
        None is returned.
    """
    return chembl_chebi.get(chembl_id)


def get_chebi_id_from_cas(cas_id):
    """Return a ChEBI ID corresponding to the given CAS ID.

    Parameters
    ----------
    cas_id : str
        The CAS ID to be converted.

    Returns
    -------
    chebi_id : str
        The ChEBI ID corresponding to the given CAS ID. If the lookup
        fails, None is returned.
    """
    return cas_chebi.get(cas_id)


def get_chebi_name_from_id(chebi_id, offline=True):
    """Return a ChEBI name corresponding to the given ChEBI ID.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID whose name is to be returned.
    offline : Optional[bool]
        If False, the ChEBI web service is invoked in case a name mapping
        could not be found in the local resource file. Default: True

    Returns
    -------
    chebi_name : str
        The name corresponding to the given ChEBI ID. If the lookup
        fails, None is returned.
    """
    chebi_id = _add_prefix(chebi_id)
    name = _obo_client.get_name_from_id(chebi_id)
    if name is None and not offline:
        return get_chebi_name_from_id_web(chebi_id)
    else:
        return name


def get_chebi_id_from_name(chebi_name):
    """Return a ChEBI ID corresponding to the given ChEBI name.

    Parameters
    ----------
    chebi_name : str
        The ChEBI name whose ID is to be returned.

    Returns
    -------
    chebi_id : str
        The ID corresponding to the given ChEBI name. If the lookup
        fails, None is returned.
    """
    return _obo_client.get_id_from_name(chebi_name)


@lru_cache(maxsize=5000)
def get_chebi_entry_from_web(chebi_id):
    """Return a ChEBI entry corresponding to a given ChEBI ID using a REST API.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID whose entry is to be returned.

    Returns
    -------
    xml.etree.ElementTree.Element
        An ElementTree element representing the ChEBI entry.
    """
    url_base = 'http://www.ebi.ac.uk/webservices/chebi/2.0/test/'
    url_fmt = url_base + 'getCompleteEntity?chebiId=%s'
    resp = requests.get(url_fmt % chebi_id)
    if resp.status_code != 200:
        logger.warning("Got bad code form CHEBI client: %s" % resp.status_code)
        return None
    tree = etree.fromstring(resp.content)
    path = 'n:Body/c:getCompleteEntityResponse/c:return'
    elem = tree.find(path, namespaces=chebi_xml_ns)
    return elem


def _get_chebi_value_from_entry(entry, key):
    if entry is None:
        return None
    path = 'c:%s' % key
    elem = entry.find(path, namespaces=chebi_xml_ns)
    if elem is not None:
        return elem.text
    return None


def get_chebi_name_from_id_web(chebi_id):
    """Return a ChEBI name corresponding to a given ChEBI ID using a REST API.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID whose name is to be returned.

    Returns
    -------
    chebi_name : str
        The name corresponding to the given ChEBI ID. If the lookup
        fails, None is returned.
    """
    entry = get_chebi_entry_from_web(chebi_id)
    return _get_chebi_value_from_entry(entry, 'chebiAsciiName')


def get_inchi_key(chebi_id):
    """Return an InChIKey corresponding to a given ChEBI ID using a REST API.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID whose InChIKey is to be returned.

    Returns
    -------
    str
        The InChIKey corresponding to the given ChEBI ID. If the lookup
        fails, None is returned.
    """
    entry = get_chebi_entry_from_web(chebi_id)
    return _get_chebi_value_from_entry(entry, 'inchiKey')


def get_primary_id(chebi_id):
    """Return the primary ID corresponding to a ChEBI ID.

    Note that if the provided ID is a primary ID, it is returned
    unchanged.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID that should be mapped to its primary equivalent.

    Returns
    -------
    str or None
        The primary ChEBI ID or None if the provided ID is neither
        primary nor a secondary ID with a primary mapping.
    """
    chebi_id = _add_prefix(chebi_id)
    if chebi_id in _obo_client.entries:
        return chebi_id
    prim_id = _obo_client.get_id_from_alt_id(chebi_id)
    return prim_id


def get_specific_id(chebi_ids):
    """Return the most specific ID in a list based on the hierarchy.

    Parameters
    ----------
    chebi_ids : list of str
        A list of ChEBI IDs some of which may be hierarchically related.

    Returns
    -------
    str
        The first ChEBI ID which is at the most specific level in the
        hierarchy with respect to the input list.
    """
    if not chebi_ids:
        return chebi_ids

    from indra.ontology.bio import bio_ontology

    def isa_cmp(a, b):
        """Compare two entries based on isa relationships for sorting."""
        if not a.startswith('CHEBI:'):
            a = 'CHEBI:%s' % a
        if not b.startswith('CHEBI:'):
            b = 'CHEBI:%s' % b
        if bio_ontology.isa('CHEBI', a, 'CHEBI', b):
            return -1
        if bio_ontology.isa('CHEBI', b, 'CHEBI', a):
            return 1
        return 0

    chebi_ids = [_add_prefix(chebi_id) for chebi_id in chebi_ids]
    chebi_id = sorted(chebi_ids, key=cmp_to_key(isa_cmp))[0]
    return chebi_id


def get_chebi_id_from_hmdb(hmdb_id):
    """Return the ChEBI ID corresponding to an HMDB ID.

    Parameters
    ----------
    hmdb_id : str
        An HMDB ID.

    Returns
    -------
    str
        The ChEBI ID that the given HMDB ID maps to or None if no mapping
        was found.
    """
    return hmdb_chebi.get(hmdb_id)


# Read resource files into module-level variables

def _read_chebi_to_pubchem():
    csv_reader = _read_resource_csv('chebi_to_pubchem.tsv')
    chebi_pubchem = {}
    pubchem_chebi = {}
    ik_matches = {}
    # Here, in case there are many possible mappings, we make it so that we
    # end up with one that has an explicit InChiKey match over one that
    # doesn't, if such a mapping is available
    for chebi_id, pc_id, ik_match in csv_reader:
        chebi_id = 'CHEBI:%s' % chebi_id
        if chebi_id not in chebi_pubchem:
            chebi_pubchem[chebi_id] = pc_id
            ik_matches[(chebi_id, pc_id)] = ik_match
        elif ik_match == 'Y' and not \
                ik_matches.get((chebi_id, chebi_pubchem[chebi_id])):
            chebi_pubchem[chebi_id] = pc_id
        if pc_id not in pubchem_chebi:
            pubchem_chebi[pc_id] = chebi_id
            ik_matches[(chebi_id, pc_id)] = ik_match
        elif ik_match == 'Y' and not \
                ik_matches.get((pubchem_chebi[pc_id], pc_id)):
            pubchem_chebi[pc_id] = chebi_id
    return chebi_pubchem, pubchem_chebi


def _read_chebi_to_chembl():
    csv_reader = _read_resource_csv('chebi_to_chembl.tsv')
    chebi_chembl = {}
    chembl_chebi = {}
    next(csv_reader)
    for row in csv_reader:
        chebi_id, chembl_id = row
        chebi_id = 'CHEBI:%s' % chebi_id
        chebi_chembl[chebi_id] = chembl_id
        chembl_chebi[chembl_id] = chebi_id
    return chebi_chembl, chembl_chebi


def _read_cas_to_chebi():
    csv_reader = _read_resource_csv('cas_to_chebi.tsv')
    cas_chebi = {}
    next(csv_reader)
    for row in csv_reader:
        cas_chebi[row[0]] = 'CHEBI:%s' % row[1]
    # These are missing from the resource but appear often, so we map
    # them manually
    extra_entries = {'24696-26-2': 'CHEBI:17761',
                     '23261-20-3': 'CHEBI:18035',
                     '165689-82-7': 'CHEBI:16618'}
    cas_chebi.update(extra_entries)
    return cas_chebi


def _read_hmdb_to_chebi():
    csv_reader = _read_resource_csv('hmdb_to_chebi.tsv')
    hmdb_chebi = {}
    next(csv_reader)
    for row in csv_reader:
        hmdb_chebi[row[0]] = 'CHEBI:%s' % row[1]
    return hmdb_chebi


def _read_resource_csv(fname):
    file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             os.pardir, 'resources', fname)
    csv_reader = read_unicode_csv(file_path, delimiter='\t')
    return csv_reader


chebi_pubchem, pubchem_chebi = _read_chebi_to_pubchem()
chebi_chembl, chembl_chebi = _read_chebi_to_chembl()
cas_chebi = _read_cas_to_chebi()
hmdb_chebi = _read_hmdb_to_chebi()
