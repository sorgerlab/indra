import indra.tools.assemble_corpus as ac
from indra.db.util import insert_pa_stmts
from indra.db.client import get_statements
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies


def make_unique_statement_set(preassembler, stmts):
    stmt_groups = preassembler._get_stmt_matching_groups(stmts)
    unique_stmts = []
    for _, duplicates in stmt_groups:
        # Get the first statement and add the evidence of all subsequent
        # Statements to it
        for stmt_ix, stmt in enumerate(duplicates):
            if stmt_ix == 0:
                first_stmt = stmt.make_generic_copy()
            first_stmt.evidence.append(stmt.uuid)
        # This should never be None or anything else
        assert isinstance(first_stmt, type(stmt))
        unique_stmts.append(first_stmt)
    return unique_stmts


def get_match_key_maps(preassembler, unique_stmts, **generate_id_map_kwargs):
    id_maps = preassembler._generate_id_maps(unique_stmts,
                                             **generate_id_map_kwargs)
    return {tuple([unique_stmts[idx].matches_key() for idx in idx_pair])
            for idx_pair in id_maps}


def process_statements(stmts, **generate_id_map_kwargs):
    stmts = ac.map_grounding(stmts)
    stmts = ac.map_sequence(stmts)
    pa = Preassembler(hierarchies)
    unique_stmts = make_unique_statement_set(pa, stmts)
    match_key_maps = get_match_key_maps(pa, unique_stmts,
                                        **generate_id_map_kwargs)
    return unique_stmts, match_key_maps


def preassemble_db_stmts(db, num_proc, *clauses):
    """Run pre-assembly on a set of statements in the database."""
    stmts = get_statements(clauses, db=db, do_stmt_count=False)
    unique_stmts, match_key_maps = process_statements(stmts, poolsize=num_proc)
    insert_pa_stmts(db, unique_stmts)
    return unique_stmts, match_key_maps


def make_graph(unique_stmts, match_key_maps):
    """Create a networkx graph of the statement and their links."""
    import networkx as nx
    g = nx.Graph()
    link_matches = {m for l in match_key_maps for m in l}
    unique_stmts_dict = {}
    for stmt in unique_stmts:
        if stmt.matches_key() in link_matches:
            g.add_node(stmt)
            unique_stmts_dict[stmt.matches_key()] = stmt

    for k1, k2 in match_key_maps:
        g.add_edge(unique_stmts_dict[k1], unique_stmts_dict[k2])

    return g
