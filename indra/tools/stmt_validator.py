import re
import json
from indra.resources import get_resource_path


class UnknownNamespace(ValueError):
    pass


class InvalidIdentifier(ValueError):
    pass


identifiers_mappings = {
    'REFSEQ_PROT': 'refseq',
    'NCBI': 'ncibgene',
    'UP': 'uniprot',
    'NONCODE': 'noncodev4.rna',
    'EGID': 'ncbigene',
    'LINCS': 'lincs.smallmolecule',
    'LNCRNADB': 'rnacentral',
    'PUBCHEM': 'pubchem.compound',
    'UPPRO': 'uniprot.chain',
    'PF': 'pfam',
    'CHEMBL': 'chembl.compound',
    'MIRBASEM': 'mirbase.mature',
    'HGNC_GROUP': 'hgnc.genefamily',
    'CTD': 'ctd.chemical',
    'IP': 'interpro',
}

non_registry = {
    'SDIS', 'SCHEM', 'SIGNOR', 'HMS-LINCS', 'NXPFA', 'GENBANK',
    'OMIM', 'LSPCI', 'ECCODE'
}

non_grounding = {
    'TEXT', 'TEXT_NORM'
}


def validate_ns(db_ns):
    identifiers_ns = identifiers_mappings.get(db_ns, db_ns.lower())
    if identifiers_ns in identifiers_registry or db_ns in non_registry \
            or db_ns in non_grounding:
        return True
    return False


def assert_valid_ns(db_ns):
    if not validate_ns(db_ns):
        raise UnknownNamespace(db_ns)


def validate_id(db_ns, db_id):
    identifiers_ns = identifiers_mappings.get(db_ns, db_ns.lower())
    if identifiers_ns in identifiers_registry:
        if re.match(identifiers_registry[identifiers_ns]['pattern'], db_id):
            return True
        else:
            return False
    elif db_ns in non_registry or db_ns in non_grounding:
        return True
    else:
        return False


def assert_valid_id(db_ns, db_id):
    if not validate_id(db_ns, db_id):
        raise InvalidIdentifier(f'{db_ns}:{db_id}')


def assert_valid_db_refs(db_refs):
    for db_ns, db_id in db_refs.items():
        assert_valid_ns(db_ns)
        assert_valid_id(db_ns, db_id)


def validate_db_refs(db_refs):
    return all(validate_ns(db_ns) and validate_id(db_ns, db_id)
               for db_ns, db_id in db_refs.items())


def load_identifiers_registry():
    with open(get_resource_path('identifiers_patterns.json'), 'r') as fh:
        return json.load(fh)


def validate_statement(stmt):
    return all(validate_db_refs(agent.db_refs)
               for agent in stmt.real_agent_list())


def assert_valid_statement(stmt):
    for agent in stmt.real_agent_list():
        assert_valid_db_refs(agent.db_refs)


identifiers_registry = load_identifiers_registry()
