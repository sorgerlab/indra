import json
from collections import defaultdict

import logging
from datetime import datetime, timedelta
from functools import wraps

import indra.tools.assemble_corpus as ac
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.db.util import insert_pa_stmts, distill_stmts_from_reading
from indra.statements import Statement
from indra.util import batch_iter

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
    def __init__(self, n_proc=1, batch_size=10000):
        self.n_proc = n_proc
        self.batch_size = batch_size
        self.pa = Preassembler(hierarchies)
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

    def _pa_batch_iter(self, db, mk_set=None):
        if mk_set is None:
            db_stmt_iter = db.select_all(db.PAStatements.json,
                                         yield_per=self.batch_size)
        else:
            db_stmt_iter = db.select_all(db.PAStatements.json,
                                         db.PAStatements.mk_hash.in_(mk_set),
                                         yield_per=self.batch_size)
        pa_stmts = (_stmt_from_json(s_json) for s_json, in db_stmt_iter)
        return batch_iter(pa_stmts, self.batch_size, return_lists=True)

    @_handle_update_table
    def create_corpus(self, db):
        # Get the statements
        _, stmt_ids = distill_stmts_from_reading(db)
        stmts = (_stmt_from_json(s_json) for s_json,
                 in db.select_all(db.RawStatements.json,
                                  db.RawStatements.uuid.in_(stmt_ids),
                                  yield_per=self.batch_size))

        # Get the set of unique statements
        hash_set = set()
        for stmt_batch in batch_iter(stmts, self.batch_size, return_lists=True):
            logger.info("Processing batch of %d statements." % len(stmt_batch))
            unique_stmts, evidence_links = \
                self._make_unique_statement_set(stmt_batch)
            new_unique_stmts = []
            for s in unique_stmts:
                s_hash = s.get_hash(shallow=True)
                if s_hash not in hash_set:
                    hash_set.add(s_hash)
                    new_unique_stmts.append(s)
            insert_pa_stmts(db, new_unique_stmts)
            db.copy('raw_unique_links', evidence_links,
                    ('pa_stmt_mk_hash', 'raw_stmt_uuid'))

        # Now get the support links between all batches.
        support_links = set()
        for i, outer_batch in enumerate(self._pa_batch_iter(db)):
            for j, inner_batch in enumerate(self._pa_batch_iter(db)):
                if i != j:
                    split_idx = len(inner_batch)
                    full_list = inner_batch + outer_batch
                    support_links |= \
                        self._get_support_links(full_list, split_idx=split_idx)
                else:
                    support_links |= self._get_support_links(inner_batch)
        db.copy('pa_support_links', support_links,
                ('supported_mk_hash', 'supporting_mk_hash'))
        return True

    @_handle_update_table
    def supplement_corpus(self, db):
        last_update = self._get_latest_updatetime(db)
        logger.info("Latest update was: %s" % last_update)

        # Get the new statements
        old_stmt_q = db.filter_query(
            db.RawStatements.uuid,
            db.RawStatements.uuid == db.RawUniqueLinks.raw_stmt_uuid
        )
        all_new_stmt_ids = (db.filter_query(db.RawStatements.uuid)
                            .except_(old_stmt_q).all())
        _, new_stmt_ids = distill_stmts_from_reading(db,
            clauses=[db.RawStatements.uuid.in_(all_new_stmt_ids)])
        new_stmts = (_stmt_from_json(s_json) for s_json,
                     in db.select_all(db.RawStatements.json,
                                      db.RawStatements.uuid.in_(new_stmt_ids),
                                      yield_per=self.batch_size))

        # Get the set of new unique statements.
        old_mk_set = {mk for mk, in db.select_all(db.PAStatements.mk_hash)}
        new_mk_set = set()
        for new_stmt_batch in batch_iter(new_stmts, self.batch_size,
                                         return_lists=True):
            logger.info("Processing batch of %d new statements."
                        % len(new_stmt_batch))
            new_unique_stmts, new_evidence_links = \
                self._make_unique_statement_set(new_stmt_batch)
            truly_new_unique_stmts = []
            for s in new_unique_stmts:
                s_hash = s.get_hash(shallow=True)
                if s_hash not in (old_mk_set | new_mk_set):
                    new_mk_set.add(s_hash)
                    truly_new_unique_stmts.append(s)
            insert_pa_stmts(db, truly_new_unique_stmts)
            db.copy('raw_unique_links', new_evidence_links,
                    ('pa_stmt_mk_hash', 'raw_stmt_uuid'))

        # Now find the new support links that need to be added.
        new_support_links = set()
        for npa_batch in self._pa_batch_iter(db, mk_set=new_mk_set):
            # Compare internally
            new_support_links |= self._get_support_links(npa_batch)

            # Compare against the other new batch statements.
            diff_new_mks = new_mk_set - {s.get_hash(shallow=True)
                                         for s in npa_batch}
            for diff_npa_batch in self._pa_batch_iter(db, mk_set=diff_new_mks):
                split_idx = len(npa_batch)
                full_list = npa_batch + diff_npa_batch
                new_support_links |= \
                    self._get_support_links(full_list, split_idx=split_idx)

            # Compare against the existing statements.
            for opa_batch in self._pa_batch_iter(db, mk_set=old_mk_set):
                split_idx = len(npa_batch)
                full_list = npa_batch + opa_batch
                new_support_links |= \
                    self._get_support_links(full_list, split_idx=split_idx)

        db.copy('pa_support_links', new_support_links,
                ('supported_mk_hash', 'supporting_mk_hash'))
        return True

    def _make_unique_statement_set(self, stmts):
        stmts = ac.map_grounding(stmts)
        stmts = ac.map_sequence(stmts)
        stmt_groups = self.pa._get_stmt_matching_groups(stmts)
        unique_stmts = []
        evidence_links = defaultdict(lambda: set())
        for _, duplicates in stmt_groups:
            # Get the first statement and add the evidence of all subsequent
            # Statements to it
            for stmt_ix, stmt in enumerate(duplicates):
                if stmt_ix == 0:
                    first_stmt = stmt.make_generic_copy()
                evidence_links[first_stmt.get_hash(shallow=True)].add(stmt.uuid)
            # This should never be None or anything else
            assert isinstance(first_stmt, type(stmt))
            unique_stmts.append(first_stmt)
        return unique_stmts, flatten_evidence_dict(evidence_links)

    def _process_statements(self, stmts, **generate_id_map_kwargs):
        unique_stmts, evidence_links = self._make_unique_statement_set(stmts)
        match_key_maps = self._get_support_links(unique_stmts,
                                                 **generate_id_map_kwargs)
        unique_stmt_dict = {stmt.get_hash(shallow=True): stmt for stmt in unique_stmts}
        return unique_stmt_dict, evidence_links, match_key_maps

    def _get_support_links(self, unique_stmts, **generate_id_map_kwargs):
        id_maps = self.pa._generate_id_maps(unique_stmts,
                                            **generate_id_map_kwargs)
        return {tuple([unique_stmts[idx].get_hash(shallow=True)
                       for idx in idx_pair])
                for idx_pair in id_maps}

    def _get_increment_links(self, unique_stmt_dict, new_unique_stmt_dict,
                             **kwargs):

        new_support_links = set()
        # Now get the list of statements the need to be compared between the
        # existing and new corpora
        old_stmt_hash_set = set(unique_stmt_dict.keys())
        new_stmt_hash_set = set(new_unique_stmt_dict.keys())
        only_old_stmts = [unique_stmt_dict[mk_hash]
                          for mk_hash in old_stmt_hash_set - new_stmt_hash_set]
        only_new_stmts = [new_unique_stmt_dict[mk_hash]
                          for mk_hash in new_stmt_hash_set - old_stmt_hash_set]
        logger.info("There were %d exclusively old statements and we are "
                    "adding %d exclusively new statements. There were %d "
                    "overlapping statements."
                    % (len(only_old_stmts), len(only_new_stmts),
                       len(new_stmt_hash_set & old_stmt_hash_set)))
        split_idx = len(only_old_stmts) + 1
        merge_stmts = only_old_stmts + only_new_stmts

        # Find the support links between the new corpora
        merge_support_links = self._get_support_links(merge_stmts,
                                                      split_idx=split_idx,
                                                      **kwargs)
        logger.info("Found %d links between old and new corpora."
                    % len(merge_support_links))

        new_support_links |= merge_support_links
        return new_support_links


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
