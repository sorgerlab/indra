from indra.statements import Phosphorylation, Agent, Evidence
from indra.db.util import NestedDict
from indra.db.util import reader_versions as rv_dict



def make_raw_statement_test_set():
    d = NestedDict()
    stmts = []
    target_sets = []
    bettered_sids = set()

    def add_stmts_to_target_set(some_stmts):
        if not target_sets:
            for stmt in some_stmts:
                target_sets.append(({stmt},
                                    {stmts.index(s) for s in some_stmts
                                     if s is not stmt}))
        else:
            old_target_sets = target_sets[:]
            target_sets.clear()
            for stmt_set, dup_set in old_target_sets:
                for stmt in some_stmts:
                    new_set = stmt_set.copy()
                    new_set.add(stmt)
                    new_dups = dup_set.copy()
                    new_dups |= {stmts.index(s) for s in some_stmts
                                 if s is not stmt}
                    target_sets.append((new_set, new_dups))
        return target_sets

    def add_content(trid, src, tcid, reader, rv_idx, rid, a, b, ev_num, copies,
                    is_target=False):
        stmts.extend(__make_test_statements(a, b, reader, ev_num, copies))
        if copies > 1:
            if ev_num is not None:
                assert stmts[-1].get_hash() == stmts[-2].get_hash()
            else:
                assert stmts[-1].get_hash() != stmts[-2].get_hash()
        rv = rv_dict[reader][rv_idx]
        r_dict = d[trid][src][tcid][reader][rv][rid]
        if ev_num is not None:
            s_hash = stmts[-1].get_hash()
            if r_dict.get(s_hash) is None:
                r_dict[s_hash] = set()
            d[trid][src][tcid][reader][rv][rid][stmts[-1].get_hash()] |= \
                {(stmts.index(s), s) for s in stmts[-copies:]}
        else:
            for s in stmts[-copies:]:
                s_hash = s.get_hash()
                if r_dict.get(s_hash) is None:
                    r_dict[s_hash] = set()
                d[trid][src][tcid][reader][rv][rid][s_hash].add(
                    (stmts.index(s), s)
                    )
        if is_target:
            global target_sets
            target_sets = add_stmts_to_target_set(stmts[-copies:])
        return

    # We produced statements a coupld of times with and old reader version
    add_content(1, 'pubmed', 1, 'reach', 0, 1, 'A1', 'B1', 1, 2)
    add_content(1, 'pubmed', 1, 'reach', 0, 1, 'A1', 'B1', 2, 1)
    add_content(1, 'pubmed', 1, 'reach', 0, 1, 'A2', 'B2', 1, 1)

    # Do it again for a new reader version.
    add_content(1, 'pubmed', 1, 'reach', 1, 2, 'A1', 'B1', 1, 2)
    add_content(1, 'pubmed', 1, 'reach', 1, 2, 'A1', 'B1', 2, 1)

    # Add some for sparser.
    add_content(1, 'pubmed', 1, 'sparser', 1, 3, 'A1', 'B1', 1, 2)
    add_content(1, 'pubmed', 1, 'sparser', 1, 3, 'A2', 'B2', 1, 1)

    # Now add statements from another source.
    add_content(1, 'pmc_oa', 2, 'reach', 0, 4, 'A1', 'B1', 1, 2)
    add_content(1, 'pmc_oa', 2, 'reach', 0, 4, 'A1', 'B1', 2, 1)
    add_content(1, 'pmc_oa', 2, 'reach', 0, 4, 'A2', 'B2', 1, 1)
    bettered_sids |= set(range(len(stmts)))

    # ...and again for a new reader version.
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A1', 'B1', 1, 2, True)
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A1', 'B1', 2, 1, True)
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A2', 'B2', 1, 1, True)
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A3', 'B3', 1, 1, True)

    # Add some results from sparser
    add_content(1, 'pmc_oa', 2, 'sparser', 1, 5, 'A1', 'B1', 1, 2, True)
    add_content(1, 'pmc_oa', 2, 'sparser', 1, 5, 'A2', 'B2', 1, 1, True)

    # Add some content for another text ref.
    add_content(2, 'pmc_oa', 3, 'sparser', 1, 6, 'A3', 'B3', 1, 1, True)
    add_content(2, 'manuscripts', 4, 'sparser', 1, 7, 'A3', 'B3', 1, 1)
    bettered_sids.add(len(stmts) - 1)

    return d, stmts, target_sets, bettered_sids


def __make_test_statements(a, b, source_api, ev_num=None, copies=1):
    stmts = []
    A = Agent(a)
    B = Agent(b)
    for i in range(copies):
        if ev_num is None:
            ev_num = i
        ev_text = "Evidence %d for %s phosphorylates %s." % (ev_num, a, b)
        ev_list = [Evidence(text=ev_text, source_api=source_api)]
        stmts.append(Phosphorylation(Agent(A), Agent(B), evidence=ev_list))
    return stmts
