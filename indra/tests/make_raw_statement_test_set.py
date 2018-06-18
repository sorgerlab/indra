from indra.statements import Phosphorylation, Agent, Evidence
from indra.db.util import NestedDict


def add_stmts_to_target_set(stmts, target_sets):
    if not target_sets:
        for stmt in stmts:
            target_sets.append(({stmt},
                                {s.uuid for s in stmts if s is not stmt}))
    else:
        old_target_sets = target_sets[:]
        target_sets = []
        for stmt_set, dup_set in old_target_sets:
            for stmt in stmts:
                new_set = stmt_set.copy()
                new_set.add(stmt)
                new_dups = dup_set.copy()
                new_dups |= {s.uuid for s in stmts if s != stmt}
                target_sets.append((new_set, new_dups))
    return target_sets


def make_raw_statement_test_set():
    d = NestedDict()
    stmts = []
    target_sets = []
    bettered_uuids = set()

    # We produced statements a coupld of times with and old reader version
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    assert stmts[-1].get_hash() == stmts[-2].get_hash()
    d[1]['pubmed'][1]['reach']['61059a-biores-e9ee36'][1][stmts[-1].get_hash()] = set(stmts[-2:])

    stmts.append(make_test_statement('A1', 'B1', 'Evidence 2 for A1 phosphorylates B1', 'reach'))
    assert stmts[-2].get_hash() != stmts[-1].get_hash()
    d[1]['pubmed'][1]['reach']['61059a-biores-e9ee36'][1][stmts[-1].get_hash()] = {stmts[-1]}

    stmts.append(make_test_statement('A2', 'B2', 'Evience 1 for A2 phosphorylates B2', 'reach'))
    d[1]['pubmed'][1]['reach']['61059a-biores-e9ee36'][1][stmts[-1].get_hash()] = {stmts[-1]}

    # Then we did it again for a new reader version.
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    assert stmts[-1].get_hash() == stmts[-2].get_hash()
    d[1]['pubmed'][1]['reach']['1.3.3-61059a-biores-'][2][stmts[-1].get_hash()] = set(stmts[-2:])
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 2 for A1 phosphorylates B1', 'reach'))
    assert stmts[-2].get_hash() != stmts[-1].get_hash()
    d[1]['pubmed'][1]['reach']['1.3.3-61059a-biores-'][2][stmts[-1].get_hash()] = {stmts[-1]}

    # Add some results from sparser
    stmts.append(make_test_statement(None, 'B1', 'Evidence 1 for A1 phosphorylates B1', 'sparser'))
    stmts.append(make_test_statement(None, 'B1', 'Evidence 1 for A1 phosphorylates B1', 'sparser'))
    assert stmts[-1].get_hash() == stmts[-2].get_hash()
    d[1]['pubmed'][1]['sparser']['sept14-linux'][3][stmts[-1].get_hash()] = set(stmts[-2:])

    stmts.append(make_test_statement(None, 'B2', 'Evidence 1 for A2 phosphoryates B2', 'sparser'))
    assert stmts[-1].get_hash() != stmts[-2].get_hash()
    d[1]['pubmed'][1]['sparser']['sept14-linux'][3][stmts[-1].get_hash()] = {stmts[-1]}

    # Now add statements from another source.
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    assert stmts[-1].get_hash() == stmts[-2].get_hash()
    d[1]['pmc_oa'][2]['reach']['61059a-biores-e9ee36'][4][stmts[-1].get_hash()] = set(stmts[-2:])

    stmts.append(make_test_statement('A1', 'B1', 'Evidence 2 for A1 phosphorylates B1', 'reach'))
    assert stmts[-2].get_hash() != stmts[-1].get_hash()
    d[1]['pmc_oa'][2]['reach']['61059a-biores-e9ee36'][4][stmts[-1].get_hash()] = {stmts[-1]}

    stmts.append(make_test_statement('A2', 'B2', 'Evience 1 for A2 phosphorylates B2', 'reach'))
    d[1]['pmc_oa'][2]['reach']['61059a-biores-e9ee36'][4][stmts[-1].get_hash()] = {stmts[-1]}

    stmts.append(make_test_statement('A2', 'B2', 'Evidence 2 for A2 phosphorylates B2', 'reach'))
    d[1]['pmc_oa'][2]['reach']['61059a-biores-e9ee36'][4][stmts[-1].get_hash()] = {stmts[-1]}

    stmts.append(make_test_statement('A3', 'B3', 'Evidence 1 for A3 phosphorylates B3', 'reach'))
    d[1]['pmc_oa'][2]['reach']['61059a-biores-e9ee36'][4][stmts[-1].get_hash()] = {stmts[-1]}
    bettered_uuids |= {s.uuid for s in stmts}

    # Then we did it again for a new reader version.
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    stmts.append(make_test_statement('A1', 'B1', 'Evidence 1 for A1 phosphorylates B1', 'reach'))
    assert stmts[-1].get_hash() == stmts[-2].get_hash()
    d[1]['pmc_oa'][2]['reach']['1.3.3-61059a-biores-'][5][stmts[-1].get_hash()] = set(stmts[-2:])
    target_sets = add_stmts_to_target_set(stmts[-2:], target_sets)

    stmts.append(make_test_statement('A1', 'B1', 'Evidence 2 for A1 phosphorylates B1', 'reach'))
    assert stmts[-2].get_hash() != stmts[-1].get_hash()
    d[1]['pmc_oa'][2]['reach']['1.3.3-61059a-biores-'][5][stmts[-1].get_hash()] = {stmts[-1]}
    target_sets = add_stmts_to_target_set(stmts[-1:], target_sets)

    stmts.append(make_test_statement('A2', 'B2', 'Evidence 2 for A2 phosphorylates B2', 'reach'))
    d[1]['pmc_oa'][2]['reach']['1.3.3-61059a-biores-'][4][stmts[-1].get_hash()] = {stmts[-1]}
    target_sets = add_stmts_to_target_set(stmts[-1:], target_sets)

    stmts.append(make_test_statement('A3', 'B3', 'Evidence 1 for A3 phosphorylates B3', 'reach'))
    d[1]['pmc_oa'][2]['reach']['1.3.3-61059a-biores-'][4][stmts[-1].get_hash()] = {stmts[-1]}
    target_sets = add_stmts_to_target_set(stmts[-1:], target_sets)

    # Add some results from sparser
    stmts.append(make_test_statement(None, 'B1', 'Evidence 1 for A1 phosphorylates B1', 'sparser'))
    stmts.append(make_test_statement(None, 'B1', 'Evidence 1 for A1 phosphorylates B1', 'sparser'))
    assert stmts[-1].get_hash() == stmts[-2].get_hash()
    d[1]['pmc_oa'][2]['sparser']['sept14-linux'][6][stmts[-1].get_hash()] = set(stmts[-2:])
    target_sets = add_stmts_to_target_set(stmts[-2:], target_sets)

    stmts.append(make_test_statement(None, 'B2', 'Evidence 1 for A2 phosphoryates B2', 'sparser'))
    assert stmts[-1].get_hash() != stmts[-2].get_hash()
    d[1]['pmc_oa'][2]['sparser']['sept14-linux'][6][stmts[-1].get_hash()] = {stmts[-1]}
    target_sets = add_stmts_to_target_set(stmts[-1:], target_sets)

    stmts.append(make_test_statement('A3', 'B3', 'Evidence 1 for A3 phosphorylates B3', 'sparser'))
    d[2]['manuscripts'][3]['sparser']['sept14-linux'][7][stmts[-1].get_hash()] = {stmts[-1]}
    target_sets = add_stmts_to_target_set(stmts[-1:], target_sets)

    stmts.append(make_test_statement('A3', 'B3', 'Evidence 1 for A3 phosphorylates B3', 'sparser'))
    d[2]['pmc_oa'][4]['sparser']['sept14-linux'][8][stmts[-1].get_hash()] = {stmts[-1]}
    bettered_uuids.add(stmts[-1].uuid)

    return d, stmts, target_sets, bettered_uuids


def make_test_statement(A, B, text, source_api):
    return Phosphorylation(Agent(A), Agent(B),
                           evidence=[Evidence(text=text, source_api=source_api)])
