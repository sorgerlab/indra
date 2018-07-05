from indra.statements import Phosphorylation, Agent, Evidence
from indra.db.util import NestedDict
from indra.db.util import reader_versions as rv_dict


def make_raw_statement_test_set():
    d = NestedDict()
    stmts = []
    target_sets = []
    bettered_sids = set()

    # Create a function which will update all possible outcome scenarios given a
    # set of some_stmts.
    def add_stmts_to_target_set(some_stmts):
        # If we don't have any target sets of statements, initialize with the
        # input statements.
        if not target_sets:
            for stmt in some_stmts:
                target_sets.append(({stmt},
                                    {stmts.index(s) for s in some_stmts
                                     if s is not stmt}))
        else:
            # Make a copy and empty the current list.
            old_target_sets = target_sets[:]
            try: # Python 3
                target_sets.clear()
            except AttributeError: # Python 2
                del target_sets[:]

            # Now for every previous scenario, pick a random possible "good"
            # statement, update the corresponding duplicate trace.
            for stmt_set, dup_set in old_target_sets:
                for stmt in some_stmts:
                    # Here we consider the possibility that each of the
                    # potential valid statements may be chosen, and record that
                    # possible alteration to the set of possible histories.
                    new_set = stmt_set.copy()
                    new_set.add(stmt)
                    new_dups = dup_set.copy()
                    new_dups |= {stmts.index(s) for s in some_stmts
                                 if s is not stmt}
                    target_sets.append((new_set, new_dups))
        return target_sets

    # Create a function to handle the creation of the metadata.
    def add_content(trid, src, tcid, reader, rv_idx, rid, a, b, ev_num, copies,
                    is_target=False):
        # Add the new statements to the over-all list.
        stmts.extend(__make_test_statements(a, b, reader, ev_num, copies))

        # If we are making multiple copies, the latest copies should have the
        # same overall hash. If it's not a copy, the hashes should be different.
        if copies > 1:
            # The above only applies if the evidence was specified to be the
            # same, otherwise it assumed the evidence, and therefore the hash,
            # is different.
            if ev_num is not None:
                assert stmts[-1].get_hash() == stmts[-2].get_hash()
            else:
                assert stmts[-1].get_hash() != stmts[-2].get_hash()

        # Populate the provenance for the dict.
        rv = rv_dict[reader][rv_idx]
        r_dict = d[trid][src][tcid][reader][rv][rid]

        # If the evidence variation was specified, the evidence in any copies is
        # identical, and they will all have the same hash. Else, the hash is
        # different and the statements need to be iterated over.
        if ev_num is not None:
            s_hash = stmts[-1].get_hash()

            # Check to see if we have a matching statment yet.
            if r_dict.get(s_hash) is None:
                r_dict[s_hash] = set()

            # Set the value
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

        # If this/these statement/s is intended to be picked up, add it/them to
        # the target sets.
        if is_target:
            global target_sets
            target_sets = add_stmts_to_target_set(stmts[-copies:])
        return

    # We produced statements a coupld of times with and old reader version
    #           trid         tcid        reader vrsn idx   distinct evidence id
    #           |  source    |  reader   |  reading id     |  number of copies
    #           |  |         |  |        |  |  Agents      |  |  Is it a target?
    add_content(1, 'pubmed', 1, 'reach', 0, 1, 'A1', 'B1', 1, 2, False)
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

    # All the statements up until now will be skipped, if all goes well.
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

    # This last statement should also be skipped, if all goes well.
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
