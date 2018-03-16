import indra.tools.assemble_corpus as ac
from indra.db.util import get_statements, insert_pa_stmts
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies


def make_unique_statement_set(stmts):
    stmt_groups = Preassembler.get_stmt_matching_groups(stmts)
    unique_stmts = []
    for _, duplicates in stmt_groups:
        # Get the first statement and add the evidence of all subsequent
        # Statements to it
        for stmt_ix, stmt in enumerate(duplicates):
            if stmt_ix == 0:
                first_stmt = stmt.get_new_copy()
            first_stmt.evidence.append(stmt.uuid)
        # This should never be None or anything else
        assert isinstance(first_stmt, type(stmt))
        unique_stmts.append(first_stmt)
    return unique_stmts


def get_match_key_maps(unique_stmts, num_procs=1):
    pass


def process_statements(stmts, num_procs=1):
    stmts = ac.map_grounding(stmts)
    stmts = ac.map_sequence(stmts)
    unique_stmts = make_unique_statement_set(stmts)
    match_key_maps = get_match_key_maps(unique_stmts, num_procs)
    return unique_stmts, match_key_maps


def preassemble_db_stmts(db, num_procs, *clauses):
    """Run pre-assembly on a set of statements in the database."""
    stmts = get_statements(clauses, db=db, do_stmt_count=False)
    pa_stmts = process_statements(stmts, num_procs)
    insert_pa_stmts(db, pa_stmts)
    return pa_stmts
