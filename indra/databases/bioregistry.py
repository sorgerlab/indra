"""This module implements a client for using namespace and identifiers
information from the Bioregistry (bioregistry.io)."""
import re
from indra.resources import load_resource_json

# Additional mappings that we need for getting a bioregistry prefix
# corresponding to an INDRA namespace, not covered by synonyms.
bioregistry_overrides = {
    'CTD': 'ctd.chemical',
    'NCBI': 'ncbi.gene',
    'NONCODE': 'noncodev4.rna',
    'PUBCHEM': 'pubchem.compound',
}


bioregistry_reverse_overrides = {v: k for k, v in bioregistry_overrides.items()}


def get_ns_from_bioregistry(bioregistry_prefix):
    """Return the INDRA namespace for the given Bioregistry prefix."""
    # If the prefix is not in Bioregistry, we return None
    if bioregistry_prefix not in registry:
        return None
    # If there is an override mapping, we return that
    mapping = bioregistry_reverse_overrides.get(bioregistry_prefix)
    if mapping:
        return mapping
    # Otherwise, we are dealing with a simple capitalization
    # conversion
    return bioregistry_prefix.upper()


def get_ns_id_from_bioregistry(bioregistry_prefix, bioregistry_id):
    """Return the INDRA namespace and ID for a Bioregistry prefix and ID."""
    # If the prefix is not in Bioregistry, we return None
    db_ns = get_ns_from_bioregistry(bioregistry_prefix)
    if not db_ns:
        return None, None
    banana = registry[bioregistry_prefix].get('banana')
    if banana:
        db_id = '%s:%s' % (banana, bioregistry_id)
    else:
        db_id = bioregistry_id
    return db_ns, db_id


def get_ns_id_from_bioregistry_curie(bioregistry_curie):
    """Return the INDRA namespace and ID for a Bioregistry CURIE."""
    # If the prefix is not in Bioregistry, we return None
    prefix, id = bioregistry_curie.split(':', maxsplit=1)
    return get_ns_id_from_bioregistry(prefix, id)


def get_bioregistry_prefix(db_ns):
    """Return the prefix for the given INDRA namespace in Bioregistry."""
    # First if there is an explicit override, we return that
    if db_ns in bioregistry_overrides:
        return bioregistry_overrides[db_ns]
    # Next, if INDRA matches a curated synonym, we return that
    if db_ns in synonym_reverse:
        return synonym_reverse[db_ns]
    # Otherwise, we check if the lowercase version of the namespace if
    # a valid prefix and return that
    if db_ns.lower() in registry:
        return db_ns.lower()
    # If none of these match, we return None
    return None


def get_bioregistry_curie(db_ns, db_id):
    """Return the Bioregistry CURIE for the given INDRA namespace and ID."""
    prefix = get_bioregistry_prefix(db_ns)
    if not prefix:
        return None
    banana = registry[prefix].get('banana')
    if banana:
        if db_id.startswith(banana):
            db_id = db_id[len(banana) + 1:]
    return '%s:%s' % (prefix, db_id)


def get_bioregistry_url(db_ns, db_id):
    """Return the Bioregistry URL for the given INDRA namespace and ID."""
    curie = get_bioregistry_curie(db_ns, db_id)
    if not curie:
        return None
    return 'https://bioregistry.io/%s' % curie


def _load_bioregistry():
    registry = load_resource_json('bioregistry.json')
    synonym_reverse = {}
    for prefix, entry in registry.items():
        # If there is a pattern we make a pre-compiled version of it
        # for faster matching.
        if 'pattern' in entry:
            pattern = entry['pattern']
            # If there is a banana, we need to add it to the pattern
            if 'banana' in entry:
                pattern = entry['banana'] + ':' + pattern
            entry['pattern_compiled'] = re.compile(pattern)
        for synonym in entry.get('synonyms', []):
            synonym_reverse[synonym] = prefix

    return registry, synonym_reverse


registry, synonym_reverse = _load_bioregistry()


