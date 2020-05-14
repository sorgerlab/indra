from indra.pipeline import AssemblyPipeline, RunnableArgument
from indra.pipeline.pipeline import jsonify_arg_input
from indra.tests.test_assemble_corpus import st1, st2, st3, st4
from indra.tools.assemble_corpus import *
from indra.preassembler.custom_preassembly import location_matches, \
    location_refinement
from indra.belief.wm_scorer import *
from indra.belief import BeliefScorer
from indra.ontology.world import world_ontology
from indra.statements import Activation


stmts = [st1, st2, st3, st4]
path_this = os.path.dirname(os.path.abspath(__file__))
test_json = os.path.join(path_this, 'pipeline_test.json')


def test_running_pipeline():
    # From json file
    ap = AssemblyPipeline.from_json_file(test_json)
    assert ap
    # AssemblyPipeline has methods for length and iteration
    assert len(ap) == 5
    for step in ap:
        assert step
    assembled_stmts = ap.run(stmts)
    assert assembled_stmts
    assert len(assembled_stmts) == 2
    # By manually adding steps
    ap2 = AssemblyPipeline()
    ap2.append(filter_no_hypothesis)
    ap2.append(map_grounding)
    ap2.append(filter_grounded_only)
    ap2.append(map_sequence)
    ap2.append(run_preassembly, return_toplevel=False)
    assembled_stmts2 = ap2.run(stmts)
    assert assembled_stmts2
    assert len(assembled_stmts2) == 2


def test_pipeline_methods():
    ap = AssemblyPipeline()
    assert len(ap) == 0
    ap.append(filter_grounded_only)
    assert len(ap) == 1
    ap.insert(0, filter_no_hypothesis)
    assert len(ap) == 2
    assert ap.steps[0] == {'function': 'filter_no_hypothesis'}
    # Append functions with arguments and runnable arguments
    ap.append(filter_by_type, Activation)
    assert len(ap) == 3
    assert ap.steps[2] == {'function': 'filter_by_type',
                           'args': [{'stmt_type': 'Activation'}]}, ap.steps[2]
    ap.append(
        run_preassembly, matches_fun=location_matches,
        refinement_fun=location_refinement, normalize_equivalences=True,
        normalize_opposites=True, normalize_ns='WM',
        belief_scoret=RunnableArgument(get_eidos_scorer),
        ontology=RunnableArgument(world_ontology))
    assert len(ap) == 4
    assert isinstance(ap.steps[3], dict)
    assert isinstance(ap.steps[3]['kwargs'], dict)
    assert len(ap.steps[3]['kwargs']) == 7
    # Run argument to get value
    assert isinstance(ap.get_argument_value({'function': 'get_eidos_scorer'}),
                      BeliefScorer)
    # Get a function object as argument
    assert ap.get_argument_value(
        {'function': 'location_matches', 'no_run': True}) == location_matches
    # Get statement type as argument
    assert ap.get_argument_value({'stmt_type': 'Activation'}) == Activation
    # Get simple argument values
    assert ap.get_argument_value('test') == 'test'
    assert ap.get_argument_value(4) == 4
    assert ap.get_argument_value(True)
    assert not ap.get_argument_value(False)
    assert ap.get_argument_value([1, 2, 3]) == [1, 2, 3]


def test_runnable_argument():
    # With args
    ra = RunnableArgument(
        get_eidos_bayesian_scorer,
        {'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]})
    assert ra
    assert ra.func_name == 'get_eidos_bayesian_scorer'
    assert ra.args == ({'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]}, )
    assert not ra.kwargs
    assert ra.to_json() == {
        'function': 'get_eidos_bayesian_scorer',
        'args': [{'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]}]}
    # With kwargs
    ra = RunnableArgument(
        get_eidos_bayesian_scorer,
        prior_counts={'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]})
    assert ra
    assert ra.func_name == 'get_eidos_bayesian_scorer'
    assert not ra.args
    assert ra.kwargs == {
        'prior_counts': {'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]}}
    assert ra.to_json() == {
        'function': 'get_eidos_bayesian_scorer',
        'kwargs': {
            'prior_counts': {
                'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]}}
    }


def test_jsonify_arg():
    assert jsonify_arg_input(4) == 4
    assert jsonify_arg_input('test') == 'test'
    assert jsonify_arg_input(True)
    assert not jsonify_arg_input(False)
    assert jsonify_arg_input([1, 2, 3]) == [1, 2, 3]
    assert jsonify_arg_input(location_matches) == {
        'function': 'location_matches', 'no_run': True}
    assert jsonify_arg_input(Activation) == {'stmt_type': 'Activation'}
    assert jsonify_arg_input(RunnableArgument(
        get_eidos_bayesian_scorer,
        {'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]})) == {
            'function': 'get_eidos_bayesian_scorer',
            'args': [{'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]}]}
