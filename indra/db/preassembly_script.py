import json
from copy import deepcopy
from collections import defaultdict

import logging
from datetime import datetime, timedelta
from functools import wraps

import indra.tools.assemble_corpus as ac
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.db.util import get_statements, insert_pa_stmts, \
    distill_stmts_from_reading
from indra.statements import stmts_from_json, Statement

logger = logging.getLogger('db_preassembly')


def _handle_update_table(func):
    @wraps(func)
    def run_and_record_update(cls, db, *args, **kwargs):
        run_datetime = datetime.utcnow()
        completed = func(cls, db, *args, **kwargs)
        if completed:
            is_corpus_init = (func.__name__ == 'create_corpus')
            db.insert('preassembly_updates', corpus_init=is_corpus_init,
                      run_datetime=run_datetime)
        return completed
    return run_and_record_update


def _stmt_from_json(stmt_json_bytes):
    return Statement._from_json(json.loads(stmt_json_bytes.decode('utf-8')))


class PreassemblyManager(object):
    """Class used to manage the preassembly pipeline"""
    def __init__(self, n_proc=1):
        self.n_proc = n_proc
        return

    def _get_latest_updatetime(self, db):
        """Get the date of the latest update."""
        update_list = db.select_all(db.PreassemblyUpdates)
        if not len(update_list):
            logger.warning("The preassembled corpus has not been initialized, "
                           "or else the updates table has not been populated.")
            return None
        return max([u.run_datetime for u in update_list])

    def _get_existing_pa_stmt_dict(self, db):
        stmt_list = db.select_all(db.PAStatements)
        return {s.mk_hash: _stmt_from_json(s.json) for s in stmt_list}

    @_handle_update_table
    def create_corpus(self, db):
        _, stmts = distill_stmts_from_reading(db, get_full_stmts=True)
        unique_stmt_dict, evidence_links, support_links = \
            process_statements(stmts, poolsize=self.n_proc)
        insert_pa_stmts(db, unique_stmt_dict.values())
        db.copy('raw_unique_links', evidence_links,
                ('pa_stmt_mk_hash', 'raw_stmt_uuid'))
        db.copy('pa_support_links', support_links,
                ('supported_mk_hash', 'supporting_mk_hash'))
        return True

    @_handle_update_table
    def supplement_corpus(self, db):
        last_update = self._get_latest_updatetime(db)
        logger.info("Latest update was: %s" % last_update)
        check_from_datetime = last_update - timedelta(hours=6)
        logger.info("Updating with statements more recent than %s."
                    % check_from_datetime)
        _, new_stmts = distill_stmts_from_reading(db, get_full_stmts=True,
                   clauses=[db.RawStatements.create_date > check_from_datetime])
        existing_stmt_dict = self._get_existing_pa_stmt_dict(db)
        new_unique_stmt_dict, new_evidence_links, new_support_links = \
            get_increment_links(existing_stmt_dict, new_stmts,
                                poolsize=self.n_proc)
        logger.info("Found %d new unique statements."
                    % len(new_unique_stmt_dict))
        only_new_hashes = set(new_unique_stmt_dict.keys()) \
                          - set(existing_stmt_dict.keys())
        insert_pa_stmts(db, [new_unique_stmt_dict[k] for k in only_new_hashes])
        db.copy('raw_unique_links', new_evidence_links,
                cols=('pa_stmt_mk_hash', 'raw_stmt_uuid'))
        db.copy('pa_support_links', new_support_links,
                cols=('supported_mk_hash', 'supporting_mk_hash'))
        return True


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
    return unique_stmts, flatten_evidence_dict(evidence_links)


def get_support_links(preassembler, unique_stmts, **generate_id_map_kwargs):
    id_maps = preassembler._generate_id_maps(unique_stmts,
                                             **generate_id_map_kwargs)
    return {tuple([unique_stmts[idx].get_shallow_hash() for idx in idx_pair])
            for idx_pair in id_maps}


def process_statements(stmts, **generate_id_map_kwargs):
    stmts = ac.map_grounding(stmts)
    stmts = ac.map_sequence(stmts)
    pa = Preassembler(hierarchies)
    unique_stmts, evidence_links = make_unique_statement_set(pa, stmts)
    match_key_maps = get_support_links(pa, unique_stmts,
                                       **generate_id_map_kwargs)
    unique_stmt_dict = {stmt.get_shallow_hash(): stmt for stmt in unique_stmts}
    return unique_stmt_dict, evidence_links, match_key_maps


def get_increment_links(unique_stmt_dict, new_stmts, **kwargs):
    # Pre-assemble the new statements.
    new_unique_stmt_dict, new_evidence_links, new_support_links = \
        process_statements(new_stmts, **kwargs)
    logger.info("Got %d new unique statments." % len(new_unique_stmt_dict))
    logger.info("Got %d new evidence links."
                % len(new_evidence_links))
    logger.info("Got %d new support links." % len(new_support_links))

    # Now get the list of statements the need to be compared between the
    # existing and new corpora

    old_stmt_hash_set = set(unique_stmt_dict.keys())
    new_stmt_hash_set = set(new_unique_stmt_dict.keys())
    only_old_stmts = [unique_stmt_dict[mk_hash]
                      for mk_hash in old_stmt_hash_set - new_stmt_hash_set]
    only_new_stmts = [new_unique_stmt_dict[mk_hash]
                      for mk_hash in new_stmt_hash_set - old_stmt_hash_set]
    logger.info("There were %d exclusively old statements and we are adding %d "
                "exclusively new statements. There were %d overlapping "
                "statements." % (len(only_old_stmts), len(only_new_stmts),
                                 len(new_stmt_hash_set & old_stmt_hash_set)))
    split_idx = len(only_old_stmts) + 1
    merge_stmts = only_old_stmts + only_new_stmts

    # Find the support links between the new corpora
    pa = Preassembler(hierarchies)
    merge_support_links = get_support_links(pa, merge_stmts,
                                            split_idx=split_idx, **kwargs)
    logger.info("Found %d links between old and new corpora."
                % len(merge_support_links))

    new_support_links |= merge_support_links
    return new_unique_stmt_dict, new_evidence_links, new_support_links


def merge_statements(unique_stmt_dict, evidence_links, support_links,
                     new_unique_stmt_dict, new_evidence_links,
                     new_support_links):
    # Update the old parameters.
    full_unique_stmt_dict = new_unique_stmt_dict.copy()
    full_unique_stmt_dict.update(unique_stmt_dict)
    logger.info("There are now %d unique statements."
                % len(full_unique_stmt_dict))

    # full_evidence_links = deepcopy(evidence_links)
    # for mk_hash, evidence_set in new_evidence_links.items():
    #     full_evidence_links[mk_hash] |= evidence_set
    full_evidence_links = evidence_links | new_evidence_links
    logger.info("There are now %d evidence links."
                % len(full_evidence_links))

    full_support_links = support_links | new_support_links
    logger.info("There are now %d support relations."
                % len(full_support_links))

    return full_unique_stmt_dict, full_evidence_links, full_support_links


def preassemble_db_stmts(db, num_proc, *clauses):
    """Run pre-assembly on a set of statements in the database."""
    stmts = get_statements(clauses, db=db, do_stmt_count=False)
    unique_stmts, evidence_links, support_links = process_statements(stmts, poolsize=num_proc)
    insert_pa_stmts(db, unique_stmts.values())
    return unique_stmts, evidence_links, support_links


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
