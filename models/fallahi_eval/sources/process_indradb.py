from indra.util import _require_python3
from indra.db.query_db_stmts import get_statements
from util import prefixed_pkl, based, basen
from indra.tools import assemble_corpus as ac


def get_db_statements():
    stmts = by_gene_role_type(agent_id='AKT1', stmt_type='Phosphorylation')
    #stmts = get_statements(do_stmt_count=False)
    ac.dump_statements(stmts, prefixed_pkl('indradb')


