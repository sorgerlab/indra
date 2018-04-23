from copy import deepcopy
from collections import defaultdict

import indra.tools.assemble_corpus as ac
from indra.db.util import insert_pa_stmts
from indra.db.client import get_statements
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies


def make_unique_statement_set(preassembler, stmts):
    stmt_groups = preassembler._get_stmt_matching_groups(stmts)
    unique_stmts = []
    evidence_links = defaultdict(lambda: set())
    for _, duplicates in stmt_groups:
        # Get the first statement and add the evidence of all subsequent
        # Statements to it
        for stmt_ix, stmt in enumerate(duplicates):
            if stmt_ix == 0:
                first_stmt = stmt.make_generic_copy()
            evidence_links[first_stmt.get_shallow_hash()].add(stmt.uuid)
        # This should never be None or anything else
        assert isinstance(first_stmt, type(stmt))
        unique_stmts.append(first_stmt)
    return unique_stmts, evidence_links


def get_match_key_maps(preassembler, unique_stmts, **generate_id_map_kwargs):
    id_maps = preassembler._generate_id_maps(unique_stmts,
                                             **generate_id_map_kwargs)
    return {tuple([unique_stmts[idx].get_shallow_hash() for idx in idx_pair])
            for idx_pair in id_maps}


def process_statements(stmts, **generate_id_map_kwargs):
    stmts = ac.map_grounding(stmts)
    stmts = ac.map_sequence(stmts)
    pa = Preassembler(hierarchies)
    unique_stmts, evidence_links = make_unique_statement_set(pa, stmts)
    match_key_maps = get_match_key_maps(pa, unique_stmts,
                                        **generate_id_map_kwargs)
    unique_stmt_dict = {stmt.get_shallow_hash(): stmt for stmt in unique_stmts}
    return unique_stmt_dict, evidence_links, match_key_maps


def merge_statements(unique_stmt_dict, evidence_links, match_key_maps,
                     new_stmts, optimize=False, **kwargs):
    # Pre-assemble the new statements.
    new_unique_stmt_dict, new_evidence_links, new_match_key_maps = \
        process_statements(new_stmts, **kwargs)

    # Now get the list of statements the need to be compared between the
    # existing and new corpora
    if optimize:
        old_stmt_hash_set = set(unique_stmt_dict.keys())
        new_stmt_hash_set = set(new_unique_stmt_dict.keys())
        only_old_stmts = [unique_stmt_dict[mk_hash]
                          for mk_hash in old_stmt_hash_set - new_stmt_hash_set]
        only_new_stmts = [new_unique_stmt_dict[mk_hash]
                          for mk_hash in new_stmt_hash_set - old_stmt_hash_set]
        split_idx = len(only_old_stmts) + 1
        merge_stmts = only_old_stmts + only_new_stmts
    else:
        raise Exception("Bad")
        split_idx = len(unique_stmt_dict) + 1
        merge_stmts = list(unique_stmt_dict.values())\
                      + list(new_unique_stmt_dict.values())

    # Find the support links between the new corpora
    pa = Preassembler(hierarchies)
    merge_match_key_maps = get_match_key_maps(pa, merge_stmts,
                                              split_idx=split_idx, **kwargs)

    # Update the old parameters.
    full_unique_stmt_dict = deepcopy(new_unique_stmt_dict)
    full_unique_stmt_dict.update(unique_stmt_dict)

    full_evidence_links = deepcopy(evidence_links)
    for mk_hash, evidence_set in new_evidence_links.items():
        evidence_links[mk_hash] |= evidence_set

    new_match_key_maps |= merge_match_key_maps
    full_match_key_maps = match_key_maps | new_match_key_maps
    return full_unique_stmt_dict, full_evidence_links, full_match_key_maps


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


def flatten_evidence_dict(ev_dict):
    return {(u_stmt_key, ev_stmt_uuid)
            for u_stmt_key, ev_stmt_uuid_set in ev_dict.items()
            for ev_stmt_uuid in ev_stmt_uuid_set}
