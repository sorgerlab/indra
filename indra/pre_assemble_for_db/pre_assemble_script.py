import indra.tools.assemble_corpus as ac
from indra.db.util import get_statements, insert_pa_stmts


def process_statements(stmts, num_procs=1):
    stmts = ac.map_grounding(stmts)
    stmts = ac.map_sequence(stmts)
    stmts = ac.run_preassembly(stmts, return_toplevel=False,
                               poolsize=num_procs)
    return stmts


def preassemble_db_stmts(db, num_procs, *clauses):
    """Run pre-assembly on a set of statements in the database."""
    stmts = get_statements(clauses, db=db, do_stmt_count=False)
    pa_stmts = process_statements(stmts, num_procs)
    insert_pa_stmts(db, pa_stmts)
    return pa_stmts
