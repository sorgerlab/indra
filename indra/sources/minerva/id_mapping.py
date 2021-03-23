from indra.databases import chebi_client
from indra.ontology.standardize import standardize_db_refs


minerva_to_indra_map = {
    'UNIPROT': 'UP',
    'NCBI_PROTEIN': 'NCBIPROTEIN',
    'REFSEQ': 'REFSEQ_PROT',
    'ENTREZ': 'EGID',
    'INTERPRO': 'IP',
    'MESH_2012': 'MESH',
    'EC': 'ECCODE',
    'PUBCHEM_SUBSTANCE': 'PUBCHEM.SUBSTANCE',
    'KEGG_COMPOUND': 'KEGG.COMPOUND'
}


def fix_id_standards(db_ns, db_id):
    if db_ns == 'CHEBI':
        if not db_id.startswith('CHEBI:'):
            db_id = f'CHEBI:{db_id}'
        db_id = chebi_client.get_primary_id(db_id)
    elif db_ns == 'HGNC' and db_id.startswith('HGNC:'):
        db_id = db_id[5:]
    return db_ns, db_id


def indra_db_refs_from_minerva_refs(refs):
    db_refs = {}
    for db_ns, db_id in refs:
        db_ns = minerva_to_indra_map[db_ns] \
            if db_ns in minerva_to_indra_map else db_ns
        db_ns, db_id = fix_id_standards(db_ns, db_id)
        db_refs[db_ns] = db_id
    # We need some special handling here for issues in the curated maps
    # If we have a specific gene grounding, remove ECCODE grounding since
    # it can incorrectly result in a family interpretation
    if 'HGNC' in db_refs:
        db_refs.pop('ECCODE', None)
    db_refs = standardize_db_refs(db_refs)
    return db_refs
