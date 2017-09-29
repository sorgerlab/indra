import json
from indra.statements import *
from indra.sources.signor import SignorProcessor, _default_csv_file
from indra.db import DatabaseManager
from indra.databases import hgnc_client

def load_stmts(db_stmt_objs):
    stmt_json_list = []
    for st_obj in db_stmt_objs:
        stmt_json_list.append(json.loads(st_obj.json.decode('utf8')))
    return stmts_from_json(stmt_json_list)

db = DatabaseManager('sqlite:///indra_test.db', 'sqlite')

#db._clear()
#sp = SignorProcessor(_default_csv_file)
#db.insert_db_stmts(sp.statements, 'signor')

def by_type(stmt_type):
    """Find all statements of a specific type."""
    db_stmts = db.select_all(db.Statements, db.Statements.type == stmt_type)
    return load_stmts(db_stmts)


def by_gene(gene_name):
    """Find all statements for a specific gene."""
    gene_id = hgnc_client.get_hgnc_id(gene_name)
    results = db.select_all([db.Statements, db.Agents],
            db.Agents.db_name == 'HGNC',
            db.Agents.db_id == gene_id,
            db.Agents.stmt_id == db.Statements.id,
    )
    db_stmts = [t[0] for t in results]
    return load_stmts(db_stmts)


def by_gene_role(gene_name, role):
    """Find all statements for a specific AGENT, ROLE.

    For example, find all regulators or targets of a specific gene.
    """
    gene_id = hgnc_client.get_hgnc_id(gene_name)
    results = db.select_all([db.Statements, db.Agents],
            db.Agents.db_name == 'HGNC',
            db.Agents.db_id == gene_id,
            db.Agents.stmt_id == db.Statements.id,
            db.Agents.role == role
    )
    db_stmts = [t[0] for t in results]
    return load_stmts(db_stmts)


def by_gene_role_type(gene_name, role, stmt_type):
    """Find stmts by gene, role, and statement type.

    For example, find phosphorylation targets of MAP2K1.
    """
    gene_id = hgnc_client.get_hgnc_id(gene_name)
    results = db.select_all([db.Statements, db.Agents],
            db.Agents.db_name == 'HGNC',
            db.Agents.db_id == gene_id,
            db.Agents.stmt_id == db.Statements.id,
            db.Agents.role == role,
            db.Statements.type == stmt_type
    )
    db_stmts = [t[0] for t in results]
    return load_stmts(db_stmts)


# Find active forms of X (filter by AGENT, TYPE)


