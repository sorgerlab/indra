from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.statements import *

ev = Evidence(source_api='bel', pmid='12345', epistemics={'direct': True},
              text='This is the evidence.')


def test_mod_condition_from():
    jd = {'mod_type': 'phosphorylation', 'residue': 'S'}
    mc = ModCondition._from_json(jd)
    assert(mc.residue == 'S')
    assert(mc.mod_type == 'phosphorylation')
    assert(mc.position is None)


def test_agent_mod_condition():
    a = Agent('MAP2K1', mods=[ModCondition('phosphorylation', 'serine', 218),
                              ModCondition('phosphorylation', 'serine', 222)])
    jd = a.to_json()
    jd2 = Agent._from_json(jd).to_json()
    assert(jd == jd2)


def test_modification():
    stmt = Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_selfmodification():
    stmt = Autophosphorylation(Agent('a'), 'Y', '1234', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_activation():
    stmt = Activation(Agent('a'), Agent('b'), 'kinase', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_amount():
    stmt = IncreaseAmount(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_active_form():
    stmt = ActiveForm(Agent('a', location='nucleus'), 'kinase', False,
                      evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_complex():
    stmt = Complex([Agent('a'), Agent('b')], evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_translocation():
    stmt = Translocation(Agent('a'), 'cytoplasm', 'nucleus', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_gap():
    stmt = Gap(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def test_gef():
    stmt = Gef(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)


def __make_support_link(supporting_stmt, supported_stmt):
    if supporting_stmt.supports is None:
        supporting_stmt.supports = []
    supporting_stmt.supports.append(supported_stmt)
    if supported_stmt.supported_by is None:
        supporting_stmt.supported_by = []
    supported_stmt.supported_by.append(supporting_stmt)
    return


def test_supports():
    stmt1 = Gap(Agent('B'), Agent('B'), evidence=[ev])
    stmt2 = Gap(Agent('a'), Agent('b'), evidence=[ev])
    __make_support_link(stmt2, stmt1)
    jd1 = stmt1.to_json()
    jd2 = stmt2.to_json()
    jds = [jd1, jd2]
    stmts = stmts_from_json(jds)
    assert(len(stmts[0].supports) == 1)
    assert(len(stmts[1].supported_by) == 1)
    assert(stmts[0].supports[0] == stmts[1])
    assert(stmts[1].supported_by[0] == stmts[0])
    jds2 = stmts_to_json(stmts)
    stmts2 = stmts_from_json(jds2)
    assert(len(stmts2[0].supports) == 1)
    assert(len(stmts2[1].supported_by) == 1)
    assert(stmts2[0].supports[0] == stmts2[1])
    assert(stmts2[1].supported_by[0] == stmts2[0])
    stmt1.to_graph()
    stmt2.to_graph()


def test_supports_missing_uuids():
    stmts = [Gap(Agent('A%d' % i), Agent('B%d' % i), evidence=[ev])
             for i in range(3)]
    for stmt in stmts[1:]:
        __make_support_link(stmts[0], stmt)
    __make_support_link(stmts[2], stmts[1])

    # Test that the example set works correctly
    output_stmts = stmts_from_json(stmts_to_json(stmts), 'error')
    assert len(output_stmts) is len(stmts),\
        "Expected %d statements, got %d." % (len(stmts), len(output_stmts))
    assert all([any([s_out.matches(s_in) for s_in in stmts])
                for s_out in output_stmts]),\
        "Output statements do not match input statements."
    all_input_supports = [s for attr in ['supports', 'supported_by']
                          for stmt in stmts for s in getattr(stmt, attr)]

    # Test the behavior when some statements are missing.
    for missing_stmt in stmts:
        print("Performing knock-out of %s." % str(missing_stmt))
        input_stmts = stmts[:]
        input_stmts.remove(missing_stmt)
        stmts_json = stmts_to_json(input_stmts)

        # Test 'handle' behavior (default)
        output_stmts = stmts_from_json(stmts_json, 'handle')
        assert len(output_stmts) == 2,\
            "Got %d stmts, not 2 when using 'handle'." % len(output_stmts)
        all_supports = [s for attr in ['supports', 'supported_by']
                        for stmt in output_stmts for s in getattr(stmt, attr)]
        unresolved_supports = [s for s in all_supports
                               if isinstance(s, Unresolved)]
        assert len(all_supports) is len(all_input_supports),\
            "Missing some supports in 'handle' behavior."
        assert 0 < len(unresolved_supports) < len(all_input_supports),\
            "No outputs are of type Unresolved for 'handle' behavior."

        # Test 'ignore' behavior
        output_stmts = stmts_from_json(stmts_json, 'ignore')
        assert len(output_stmts) == 2,\
            "Got %d stmts, not 2 when using 'ignore'." % len(output_stmts)
        all_supports = [s for attr in ['supports', 'supported_by']
                        for stmt in output_stmts for s in getattr(stmt, attr)]
        exp_num_remaining = len(all_input_supports) - len(unresolved_supports)
        assert len(all_supports) == exp_num_remaining,\
            ("Expected %d remaining support stmts, but got %d using 'ignore'."
             % (exp_num_remaining, len(all_supports)))

        # Test 'error' behavior
        try:
            output_stmts = stmts_from_json(stmts_json, 'error')
            assert False, "Failed to error when passing partial set of stmts."
        except UnresolvedUuidError:
            pass
    return
