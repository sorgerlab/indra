import json
from indra.statements import *
from indra.sources.signor import SignorProcessor, _default_csv_file
from indra.db import DatabaseManager, get_primary_db
from indra.databases import hgnc_client


db = get_primary_db()

#db = DatabaseManager('sqlite:///indra_test.db', 'sqlite')
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
    stmts = get_statements(*clauses)
    #results = db.select_all([db.Statements, db.Agents], *clauses)
    #db_stmts = [t[0] for t in results]
    return stmts


def get_statements(*clauses, count=1000):
    stmts = []
    q = db.filter_query('statements', *clauses)
    print("Counting statements...")
    num_stmts = q.count()
    print("Total of %d statements" % num_stmts)
    db_stmts = q.yield_per(count)
    subset = []
    total_counter = 0
    for stmt in db_stmts:
        subset.append(stmt)
        if len(subset) == count:
            stmts.extend(stmts_from_db_list(subset))
            subset = []
        total_counter += 1
        if total_counter % count == 0:
            print("%d of %d" % (total_counter, num_stmts))
    stmts.extend(stmts_from_db_list(subset))
    return stmts


def stmts_from_db_list(db_stmt_objs):
    stmt_json_list = []
    for st_obj in db_stmt_objs:
        stmt_json_list.append(json.loads(st_obj.json.decode('utf8')))
    return stmts_from_json(stmt_json_list)

