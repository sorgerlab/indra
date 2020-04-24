from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import datetime
import jsonschema
from indra.statements import *
from .test_json_schema import schema

ev = Evidence(source_api='bel', pmid='12345', epistemics={'direct': True},
              text='This is the evidence.')


def test_mod_condition_from():
    jd = {'mod_type': 'phosphorylation', 'residue': 'S'}
    mc = ModCondition._from_json(jd)
    assert mc.residue == 'S'
    assert mc.mod_type == 'phosphorylation'
    assert mc.position is None


def test_agent_mod_condition():
    a = Agent('MAP2K1', mods=[ModCondition('phosphorylation', 'serine', 218),
                              ModCondition('phosphorylation', 'serine', 222)])
    jd = a.to_json()
    jd2 = Agent._from_json(jd).to_json()
    assert jd == jd2


def test_modification():
    stmt = Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_selfmodification():
    stmt = Autophosphorylation(Agent('a'), 'Y', '1234', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_activation():
    stmt = Activation(Agent('a'), Agent('b'), 'kinase', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_amount():
    stmt = IncreaseAmount(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_active_form():
    stmt = ActiveForm(Agent('a', location='nucleus'), 'kinase', False,
                      evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_complex():
    stmt = Complex([Agent('a'), Agent('b')], evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_translocation():
    stmt = Translocation(Agent('a'), 'cytoplasm', 'nucleus', evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_gap():
    stmt = Gap(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_gef():
    stmt = Gef(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    stmt.to_graph()
    jd2 = Statement._from_json(jd).to_json()
    assert jd == jd2


def test_influence():
    ev1 = Event(Concept('inorganic fertilizer'),
                delta=QualitativeDelta(polarity=1, adjectives=['serious']),
                context=WorldContext(geo_location=RefContext('x')))
    ev2 = Event(Concept('farm sizes'),
                delta=QualitativeDelta(polarity=1, adjectives=['significant']))
    stmt = Influence(ev1, ev2)
    jd = stmt.to_json()
    assert 'sbo' not in jd['subj']
    assert 'sbo' not in jd['obj']
    stmt.to_graph()
    st_deserialize = Statement._from_json(jd)
    assert st_deserialize.subj.delta.polarity == 1
    assert st_deserialize.obj.delta.adjectives == ['significant']
    assert st_deserialize.subj.context.geo_location.name == 'x', \
        st_deserialize.subj.context
    jd2 = st_deserialize.to_json()
    assert jd == jd2, (jd, jd2)
    jsonschema.validate([json.loads(json.dumps(jd))], schema)


def test_association():
    ev1 = Event(Concept('food'), 
                delta=QualitativeDelta(polarity=-1, adjectives=['insufficient']))
    ev2 = Event(Concept('hunger'),
                delta=QualitativeDelta(polarity=1, adjectives=None))
    stmt = Association([ev1, ev2])
    jd = stmt.to_json()
    stmt.to_graph()
    assert 'members' in jd
    st_deserialize = Statement._from_json(jd)
    assert isinstance(st_deserialize.members[0], Event)
    assert st_deserialize.members[0].delta.polarity == -1
    assert st_deserialize.members[0].delta.adjectives == ['insufficient']
    jd2 = st_deserialize.to_json()
    assert jd == jd2


def test_migration():
    m = Migration(Concept('migration'), 
                  QuantitativeState('person', 500, 'absolute'),
                  MovementContext(locations=[
                      {'location': RefContext('South Sudan'), 'role': 'origin'},
                      {'location': RefContext('Ethiopia'), 'role': 'destination'}],
                      time=TimeContext(text='beginning of 2019')))
    jd = m.to_json()
    m.to_graph()
    assert jd
    st_deserialize = Statement._from_json(jd)
    assert isinstance(st_deserialize, Migration)
    assert isinstance(st_deserialize.delta, QuantitativeState)
    assert isinstance(st_deserialize.context, MovementContext)
    jd2 = st_deserialize.to_json()
    assert jd == jd2


def __make_support_link(supporting_stmt, supported_stmt):
    supporting_stmt.supports.append(supported_stmt)
    supported_stmt.supported_by.append(supporting_stmt)
    return


def test_supports():
    stmt1 = Gap(Agent('B'), Agent('B'), evidence=[ev])
    stmt2 = Gap(Agent('a'), Agent('b'), evidence=[ev])
    __make_support_link(stmt1, stmt2)
    jd1 = stmt1.to_json()
    jd2 = stmt2.to_json()
    jds = [jd1, jd2]
    stmts = stmts_from_json(jds)
    assert len(stmts[0].supports) == 1
    assert len(stmts[1].supported_by) == 1
    assert stmts[0].supports[0] == stmts[1]
    assert stmts[1].supported_by[0] == stmts[0]
    jds2 = stmts_to_json(stmts)
    stmts2 = stmts_from_json(jds2)
    assert len(stmts2[0].supports) == 1
    assert len(stmts2[1].supported_by) == 1
    assert stmts2[0].supports[0] == stmts2[1]
    assert stmts2[1].supported_by[0] == stmts2[0]
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
    all_input_supports = [s for stmt in stmts
                          for s in (stmt.supports + stmt.supported_by)]
    print("Total input supports: %d." % len(all_input_supports))

    # Test the behavior when some statements are missing.
    for missing_stmt in stmts:
        print("Performing knock-out of %s." % str(missing_stmt))
        input_stmts = stmts[:]
        input_stmts.remove(missing_stmt)
        exp_num_supports_removed = (len(missing_stmt.supports)
                                    + len(missing_stmt.supported_by))
        print("Knockout removed %d supports." % exp_num_supports_removed)
        stmts_json = stmts_to_json(input_stmts)

        # Test 'handle' behavior (default)
        output_stmts = stmts_from_json(stmts_json, 'handle')
        assert len(output_stmts) == 2,\
            "Got %d stmts, not 2 when using 'handle'." % len(output_stmts)
        all_supports = [s for stmt in output_stmts
                        for s in (stmt.supports + stmt.supported_by)]
        print("Number of supports for 'handle': %d." % len(all_supports))
        unresolved_supports = [s for s in all_supports
                               if isinstance(s, Unresolved)]
        unresolved_uuids = get_unresolved_support_uuids(output_stmts)
        assert unresolved_uuids == {s.uuid for s in unresolved_supports}
        print("Number of unresolved supports: %d." % len(unresolved_supports))
        exp_num_handle_supports = (len(all_input_supports)
                                   - exp_num_supports_removed)
        print("Expected number of supports: %d." % exp_num_handle_supports)
        assert len(all_supports) is exp_num_handle_supports,\
            ("Expected %d supports in 'handle' behavior, got %d."
             % (exp_num_handle_supports, len(all_supports)))
        assert len(unresolved_supports) is exp_num_supports_removed,\
            ("Expected %d output supports of type Unresolved "
             "for 'handle' behavior, but got %d."
             % (exp_num_supports_removed, len(unresolved_supports)))

        # Test 'ignore' behavior
        output_stmts = stmts_from_json(stmts_json, 'ignore')
        assert len(output_stmts) == 2,\
            "Got %d stmts, not 2 when using 'ignore'." % len(output_stmts)
        all_supports = [s for stmt in output_stmts
                        for s in (stmt.supports + stmt.supported_by)]
        print("Number of supports for 'ignore': %d." % len(all_supports))
        exp_num_ignore_supports = (exp_num_handle_supports
                                   - len(unresolved_supports))
        assert len(all_supports) == exp_num_ignore_supports,\
            ("Expected %d remaining support stmts, but got %d using 'ignore'."
             % (exp_num_ignore_supports, len(all_supports)))

        # Test 'error' behavior
        try:
            output_stmts = stmts_from_json(stmts_json, 'error')
            assert False, "Failed to error when passing partial set of stmts."
        except UnresolvedUuidError:
            pass
    return


def test_belief():
    stmt = Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev])
    jd = stmt.to_json()
    assert 'belief' in jd
    assert jd['belief'] == 1
    stmt.belief  = 0.5
    jd2 = Statement._from_json(stmt.to_json()).to_json()
    assert jd2['belief'] == 0.5


def test_hash():
    # Check that the original matches hash is preserved to and from JSON
    stmt = Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev])
    sh = stmt.get_hash()
    jd = stmt.to_json()
    stmt_ret = Statement._from_json(jd)
    assert stmt_ret._shallow_hash == sh

    # Check that the JSON matches hash is set but when refreshed, is replaced
    jd['matches_hash'] = '1234'
    stmt_ret = Statement._from_json(jd)
    assert stmt_ret._shallow_hash == 1234, stmt_ret._shallow_hash
    assert stmt_ret.get_hash() == 1234
    assert stmt_ret.get_hash(refresh=True) == sh

    # Check handling of invalid matches hash in JSON
    jd['matches_hash'] = 'abcde'
    stmt_ret = Statement._from_json(jd)
    assert stmt_ret._shallow_hash is None
    assert stmt_ret.get_hash() == sh


def test_time_context():
    tc = TimeContext(text='2018',
                     start=datetime.datetime(2018, 1, 1, 0, 0),
                     end=datetime.datetime(2019, 1, 1, 0, 0),
                     duration=(365 * 86400))
    jd = tc.to_json()
    assert jd['text'] == '2018'
    assert jd['start'] == '2018-01-01T00:00'
    assert jd['end'] == '2019-01-01T00:00'
    assert jd['duration'] == 365 * 86400

    assert TimeContext.from_json(jd).__dict__ == tc.__dict__


def test_ref_context():
    rc1 = RefContext(name='x', db_refs={'y': '1', 'z': '2'})
    rc2 = RefContext(name='x')
    rc3 = RefContext(db_refs={'y'})
    rj1 = rc1.to_json()
    assert rj1['name'] == 'x'
    assert rj1['db_refs'] == {'y': '1', 'z': '2'}
    assert rc1.to_json() == RefContext.from_json(rc1.to_json()).to_json()
    assert rc2.to_json() == RefContext.from_json(rc2.to_json()).to_json()
    assert rc3.to_json() == RefContext.from_json(rc3.to_json()).to_json()


def test_bio_context():
    rc1 = RefContext(name='x', db_refs={'y': '1', 'z': '2'})
    rc2 = RefContext(name='x')
    rc3 = RefContext(db_refs={'y'})
    bc = BioContext(location=rc1, cell_line=rc2, species=rc3)
    bcj = bc.to_json()
    assert bcj['type'] == 'bio'
    assert bcj['location']['name'] == 'x'
    assert bcj['location']['db_refs'] == {'y': '1', 'z': '2'}
    assert bc.to_json() == Context.from_json(bc.to_json()).to_json()


def test_world_context():
    gl = RefContext(name='x', db_refs={'y': '1', 'z': '2'})
    tc = TimeContext(text='2018')
    wc = WorldContext(time=tc, geo_location=gl)
    wcj = wc.to_json()
    assert wcj['type'] == 'world'
    assert wcj['time']['text'] == '2018'
    assert wcj['geo_location']['name'] == 'x'
    assert wcj['geo_location']['db_refs'] == {'y': '1', 'z': '2'}
    assert wc.to_json() == Context.from_json(wc.to_json()).to_json()


def test_movement_context():
    origin = RefContext(name='x', db_refs={'a': '1', 'b': '2'})
    interim = RefContext(name='y', db_refs={'c': '3', 'd': '4'})
    destination = RefContext(name='z', db_refs={'e': '5', 'f': '6'})
    tc = TimeContext(text='2018')
    mc = MovementContext(
        locations=[
            {'location': origin, 'role': 'origin'},
            {'location': interim, 'role': 'interim'},
            {'location': destination, 'role': 'destination'}],
        time=tc)
    mcj = mc.to_json()
    assert mcj['type'] == 'movement'
    assert mcj['locations'][0]['location']['name'] == 'x'
    assert mcj['locations'][0]['location']['db_refs'] == {'a': '1', 'b': '2'}
    assert mcj['locations'][0]['role'] == 'origin'
    assert mcj['time']['text'] == '2018'
    assert mc.to_json() == Context.from_json(mcj).to_json()


def test_evidence_context():
    gl = RefContext(name='x', db_refs={'y': '1', 'z': '2'})
    tc = TimeContext(text='2018')
    wc = WorldContext(time=tc, geo_location=gl)
    ev = Evidence(pmid='1', text='x', annotations={'a': '2'},
                  context=wc)
    evj = ev.to_json()
    assert evj['context']['type'] == 'world'
    assert evj['context']['geo_location']['name'] == 'x'
    assert evj['pmid'] == '1'
    assert evj['annotations'] == {'a': '2'}
    assert ev.to_json() == Evidence._from_json(ev.to_json()).to_json()


def test_file_serialization():
    stmt = IncreaseAmount(Agent('a'), Agent('b'), evidence=[ev])
    stmts_to_json_file([stmt], 'test_indra_stmts.json')
    stmts = stmts_from_json_file('test_indra_stmts.json')
    assert stmts[0].matches(stmt)


def test_file_serialization_json_lines():
    stmt = IncreaseAmount(Agent('a'), Agent('b'), evidence=[ev])
    stmts_to_json_file([stmt], 'test_indra_stmts.json', format='jsonl')
    stmts = stmts_from_json_file('test_indra_stmts.json', format='jsonl')
    assert stmts[0].matches(stmt)
