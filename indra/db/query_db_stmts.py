import json
from indra.statements import *
from indra.sources.signor import SignorProcessor, _default_csv_file
from indra.db import DatabaseManager
from indra.databases import hgnc_client


db = DatabaseManager('sqlite:///indra_test.db', 'sqlite')


#db._clear()
#sp = SignorProcessor(_default_csv_file)
#db.insert_db_stmts(sp.statements, 'signor')


def by_gene_role_type(gene_name=None, role=None, stmt_type=None):
    if not (gene_name or role or stmt_type):
        raise ValueError('At least one of gene_name, role, or stmt_type '
                         'must be specified.')
    clauses = []
    if gene_name:
        gene_id = hgnc_client.get_hgnc_id(gene_name)
        clauses.extend([db.Agents.db_name == 'HGNC',
                        db.Agents.db_id == gene_id])
    if role:
        clauses.append(db.Agents.role == role)
    if gene_name or role:
        clauses.append(db.Agents.stmt_id == db.Statements.id)
    if stmt_type:
        clauses.append(db.Statements.type == stmt_type)
    results = db.select_all([db.Statements, db.Agents], *clauses)
    db_stmts = [t[0] for t in results]
    return load_stmts(db_stmts)


def load_stmts(db_stmt_objs):
    stmt_json_list = []
    for st_obj in db_stmt_objs:
        stmt_json_list.append(json.loads(st_obj.json.decode('utf8')))
    return stmts_from_json(stmt_json_list)

