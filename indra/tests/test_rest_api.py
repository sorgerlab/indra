import json
from datetime import datetime
from copy import deepcopy
from nose.plugins.attrib import attr
from os import path
from rest_api.api import api
from indra.statements import *


HERE = path.dirname(path.abspath(__file__))

a = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'})
b = Agent('b', db_refs={'UP': 'P15056', 'TEXT': 'b'})
c = Agent('c', db_refs={'FPLX': 'XXX', 'TEXT': 'c'})
d = Agent('d', db_refs={'TEXT': 'd'})
e = Agent('e', db_refs={'CHEBI': 'CHEBI:1234', 'TEXT': 'e'})
f = Agent('b', db_refs={'UP': 'P28028', 'TEXT': 'b'})
g = Agent('g', db_refs={'FPLX': 'ERK'})
h = Agent('g', mods=[ModCondition('phosphorylation', 'S', '202')],
          mutations=[MutCondition('858', 'L', 'R')],
          activity=ActivityCondition('activity', True), location='nucleus',
          bound_conditions=[BoundCondition(a), BoundCondition(b)])
i = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'},
          bound_conditions=[BoundCondition(d)])
j = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'},
          bound_conditions=[BoundCondition(b)])
k = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'},
          bound_conditions=[BoundCondition(f)])
l = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'},
          bound_conditions=[BoundCondition(a)])

mapk1 = Agent('MAPK1', db_refs={'HGNC':'6871', 'UP':'P28482'})
erk = Agent('ERK', db_refs={'FPLX': 'ERK'})
st1 = Phosphorylation(
    a, b, evidence=[Evidence(text='a->b', source_api='assertion')])
st2 = Phosphorylation(
    a, d, evidence=[Evidence(text='a->d', source_api='assertion')])
st3 = Phosphorylation(
    c, d, evidence=[Evidence(text='c->d', source_api='assertion')])
st4 = Phosphorylation(
    b, e, evidence=[Evidence(text='b->e', source_api='assertion')])
st5 = Phosphorylation(
    None, b, evidence=[Evidence(text='->b', source_api='assertion')])
st6 = Phosphorylation(
    None, d, evidence=[Evidence(text='->d', source_api='assertion')])
st7 = Phosphorylation(
    None, e, evidence=[Evidence(text='->e', source_api='assertion')])
st8 = Phosphorylation(
    b, f, evidence=[Evidence(text='b->f', source_api='assertion')])
st9 = Phosphorylation(
    None, f, evidence=[Evidence(text='->f', source_api='assertion')])
st10 = Phosphorylation(
    None, g, evidence=[Evidence(text='->g', source_api='assertion')])
st11 = Phosphorylation(
    None, h, evidence=[Evidence(text='->h', source_api='assertion')])
st12 = Phosphorylation(
    a, b, evidence=[Evidence(epistemics={'direct': True})])
st13 = Phosphorylation(
    a, b, evidence=[Evidence(epistemics={'direct': False})])
st14 = Activation(a, b, 'activity')
st15 = Activation(a, b, 'kinase')
st14.supports = [st15]
st15.supported_by = [st14]
st16 = Phosphorylation(a, mapk1)
st17 = Phosphorylation(a, erk)
st18 = Phosphorylation(a, i)
st19 = Phosphorylation(a, j)
st20 = Phosphorylation(a, k)
st21 = Phosphorylation(a, l)
st1.belief = 0.9
st2.belief = 0.8
st3.belief = 0.7


def _call_api(method, route, *args, **kwargs):
    tc = api.app.test_client()
    req_meth = getattr(tc, method)
    start = datetime.now()
    print("Submitting request to '%s' at %s." % ('/' + route, start))
    print("\targs:", args)
    print("\tkwargs:", kwargs)
    res = req_meth(route, *args, **kwargs)
    end = datetime.now()
    print("Got result with %s at %s after %s seconds."
          % (res.status_code, end, (end-start).total_seconds()))
    assert res.status_code == 200, res.status_code
    return res


def _post_stmts_preassembly(stmts, route, **kwargs):
    stmts_json = stmts_to_json(stmts)
    req_json = {'statements': stmts_json}
    if kwargs:
        req_json.update(kwargs)
    res = _call_api('post', route, json=req_json)
    res_json = json.loads(res.get_data())
    assert 'statements' in res_json
    st_out = stmts_from_json(res_json.get('statements'))
    return st_out


def test_filter_grounded_only():
    route = 'preassembly/filter_grounded_only'
    st_out = _post_stmts_preassembly([st1, st4], route)
    assert len(st_out) == 2
    st_out = _post_stmts_preassembly([st3], route)
    assert len(st_out) == 0

    # Do we filter out a statement with an ungrounded bound condition?
    st_out = _post_stmts_preassembly([st18], route)
    assert len(st_out) == 0

    # When we request to remove ungrounded bound conditions, do we?
    st18_copy = deepcopy(st18)
    assert len(st18_copy.sub.bound_conditions) == 1
    st_out = _post_stmts_preassembly([st18_copy], route, remove_bound=True)
    assert len(st_out[0].sub.bound_conditions) == 0

    # When we request to remove ungrounded bound conditions, do we leave
    # grounded bound conditions in place?
    st19_copy = deepcopy(st19)
    assert len(st19_copy.sub.bound_conditions) == 1
    st_out = _post_stmts_preassembly([st19_copy], route, remove_bound=True)
    assert len(st_out[0].sub.bound_conditions) == 1

    # Do we filter out a statement with an grounded bound condition?
    st_out = _post_stmts_preassembly([st19], route)
    assert len(st_out) == 1


def test_filter_grounded_only_score():
    route = 'preassembly/filter_grounded_only'
    c1 = Event(Concept('x', db_refs={'a': [('x', 0.5), ('y', 0.8)]}))
    c2 = Event(Concept('x', db_refs={'a': [('x', 0.7), ('y', 0.9)]}))
    st1 = Influence(c1, c2)
    st_out = _post_stmts_preassembly([st1], route)
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly([st1], route, score_threshold=0.4)
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly([st1], route, score_threshold=0.6)
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly([st1], route, score_threshold=0.85)
    assert len(st_out) == 0
    st_out = _post_stmts_preassembly([st1], route, score_threshold=0.95)
    assert len(st_out) == 0
    c3 = Event(Concept('x', db_refs={'a': []}))
    st2 = Influence(c1, c3)
    assert len(_post_stmts_preassembly([st2], route)) == 0


def test_filter_uuid_list():
    route = 'preassembly/filter_uuid_list'
    st_out = _post_stmts_preassembly([st1, st4], route, uuids=[st1.uuid])
    assert len(st_out) == 1


def test_filter_genes_only():
    route = 'preassembly/filter_genes_only'
    st_out = _post_stmts_preassembly([st1, st5], route)
    assert len(st_out) == 2
    st_out = _post_stmts_preassembly([st6, st7], route)
    assert len(st_out) == 0
    st_out = _post_stmts_preassembly([st4], route)
    assert len(st_out) == 0
    st_out = _post_stmts_preassembly([st3], route, specific_only=True)
    assert len(st_out) == 0

    # Can we remove statements with non-gene bound conditions?
    st_out = _post_stmts_preassembly([st18], route)
    assert len(st_out) == 0
    st_out = _post_stmts_preassembly([st18], route, remove_bound=False)
    assert len(st_out) == 0

    # Can we remove non-gene bound conditions?
    st18_copy = deepcopy(st18)
    assert len(st18_copy.sub.bound_conditions) == 1
    st_out = _post_stmts_preassembly([st18_copy], route, remove_bound=True)
    assert len(st_out[0].sub.bound_conditions) == 0


def test_filter_human_only():
    route = 'preassembly/filter_human_only'
    st_out = _post_stmts_preassembly([st1, st5], route)
    assert len(st_out) == 2
    st_out = _post_stmts_preassembly([st8, st9], route)
    assert len(st_out) == 0

    # Can we filter out statements with bound conditions grounded to non-human
    # genes?
    st_out = _post_stmts_preassembly([st20], route, remove_bound=False)
    assert len(st_out) == 0

    # When we do such filtering, do we keep statements bounded to human genes?
    st_out = _post_stmts_preassembly([st21], route, remove_bound=False)
    assert len(st_out) == 1

    # Can we remove bound conditions grounded to non-human genes?
    st_out = _post_stmts_preassembly([st20], route, remove_bound=True)
    assert len(st_out) == 1
    assert len(st_out[0].sub.bound_conditions) == 0

    # When we do so, do we keep bound conditions not grounded to non-human
    # genes?
    st_out = _post_stmts_preassembly([st21], route, remove_bound=False)
    assert len(st_out) == 1
    assert len(st_out[0].sub.bound_conditions) == 1


def test_filter_gene_list_one():
    route = 'preassembly/filter_gene_list'
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['a'], policy='one')
    assert len(st_out) == 2
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['a'], policy='all')
    assert len(st_out) == 0
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['a', 'b'], policy='all')
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['a', 'b'], policy='invalid')
    assert len(st_out) == 2

    # Can we exclude a statement with a bound condition agent not on the filter
    # list?
    st_out = _post_stmts_preassembly(
        [st18], route, gene_list=['a', 'b', 'd'], policy='all')
    # All genes in the list
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly(
        [st18], route, gene_list=['a', 'b'], policy='all')
    # Bound condition for sub not in list
    assert len(st_out) == 0
    st_out = _post_stmts_preassembly(
        [st18], route, gene_list=['a', 'b'], policy='one')
    # Bound condition for sub not in list but we only need to match one
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly(
        [st18], route, gene_list=['d'], policy='one')
    # Only the bound condition is in filter list
    assert len(st_out) == 1

    # Can we remove bound conditions that are not in the filter list?
    st_out = _post_stmts_preassembly(
        [st18], route, gene_list=['a', 'b', 'd'], policy='all',
        remove_bound=True)
    assert len(st_out[0].sub.bound_conditions) == 1
    st_out = _post_stmts_preassembly(
        [st18], route, gene_list=['a', 'b'], policy='all', remove_bound=True)
    assert len(st_out[0].sub.bound_conditions) == 0


def test_filter_gene_list_invert():
    route = 'preassembly/filter_gene_list'
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['a'], policy='one', invert=True)
    assert len(st_out) == 0
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['d'], policy='one', invert=True)
    assert len(st_out) == 1
    assert st_out[0].sub.name == 'b'
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['a', 'd'], policy='all', invert=True)
    assert len(st_out) == 1
    assert st_out[0].sub.name == 'b'
    st_out = _post_stmts_preassembly(
        [st1, st2], route, gene_list=['a', 'b', 'd'], policy='all',
        invert=True)
    assert len(st_out) == 0


def test_filter_gene_list_families():
    route = 'preassembly/filter_gene_list'
    st_out = _post_stmts_preassembly(
        [st16, st17], route, gene_list=['MAPK1'], policy='one',
        allow_families=False)
    assert len(st_out) == 1
    assert st_out[0].equals(st16)
    st_out = _post_stmts_preassembly(
        [st16, st17], route, gene_list=['MAPK1'], policy='one',
        allow_families=True)
    assert len(st_out) == 2
    assert any([st16.equals(st) for st in st_out])
    assert any([st17.equals(st) for st in st_out])


def test_run_preassembly():
    route = 'preassembly/run_preassembly'
    st_out = _post_stmts_preassembly([st1, st3, st5, st6], route)
    assert len(st_out) == 2


def test_run_preassembly_all_stmts():
    route = 'preassembly/run_preassembly'
    st_out = _post_stmts_preassembly(
        [st1, st3, st5, st6], route, return_toplevel=False)
    assert len(st_out) == 4


def test_run_preassembly_concepts():
    route = 'preassembly/run_preassembly'
    rainfall = Event(Concept('rain', db_refs={
        'WM': 'wm/concept/indicator_and_reported_property/weather/rainfall'}))
    flooding_1 = Event(Concept('flood', db_refs={
        'WM': 'wm/x/y/z/flooding'}))
    flooding_2 = Event(Concept('flooding', db_refs={
        'WM': 'wm/a/b/c/flooding'}))
    st_out = _post_stmts_preassembly(
        [Influence(rainfall, flooding_1), Influence(rainfall, flooding_2)],
        route, normalize_ns='WM', normalize_equivalences=True,
        belief_scorer='wm', ontology='wm')
    assert len(st_out) == 2, st_out


def test_preassembly_wm_scorer():
    route = 'preassembly/run_preassembly'
    ev = Evidence(source_api='eidos',
                  annotations={'found_by': 'dueToSyntax2-Causal'})
    st = Influence(Event(Concept('x')), Event(Concept('y')), evidence=[ev])
    stmts_json = stmts_to_json([st])
    st_out = _post_stmts_preassembly([st], route, belief_scorer='wm')
    stmt = st_out[0]
    assert stmt.belief == 0.8142, stmt.belief


def test_expand_families():
    route = 'preassembly/expand_families'
    st_out = _post_stmts_preassembly([st10], route)
    assert len(st_out) == 2


def test_strip_agent_context():
    route = 'preassembly/strip_agent_context'
    st_out = _post_stmts_preassembly([st11], route)
    assert len(st_out) == 1
    assert not st_out[0].sub.mods
    assert not st_out[0].sub.mutations
    assert not st_out[0].sub.bound_conditions
    assert not st_out[0].sub.activity
    assert not st_out[0].sub.location


def test_filter_direct():
    route = 'preassembly/filter_direct'
    st_out = _post_stmts_preassembly([st12], route)
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly([st13], route)
    assert len(st_out) == 0


def test_filter_belief():
    route = 'preassembly/filter_belief'
    st_out = _post_stmts_preassembly(
        [st1, st2, st3], route, belief_cutoff=0.75)
    assert len(st_out) == 2


def test_reduce_activities():
    route = 'preassembly/reduce_activities'
    st_out = _post_stmts_preassembly([st14, st15], route)
    assert st_out[0].obj_activity == 'kinase'
    assert st_out[1].obj_activity == 'kinase'


def test_filter_source():
    route = 'preassembly/filter_evidence_source'
    ev1 = Evidence(source_api='bel')
    ev2 = Evidence(source_api='biopax')
    ev3 = Evidence(source_api='reach')
    st1 = Activation(Agent('a'), Agent('b'), evidence=[ev3])
    st2 = Activation(Agent('a'), Agent('b'), evidence=[ev1, ev2])
    st3 = Activation(Agent('a'), Agent('b'), evidence=[ev1, ev3])
    st_out = _post_stmts_preassembly([st1, st2], route, source_apis=['reach'],
                                     policy='one')
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly([st1, st2, st3], route,
                                     source_apis=['reach'], policy='all')
    assert (len(st_out) == 2)
    st_out = _post_stmts_preassembly(
        [st1, st2, st3], route, source_apis=['bel', 'biopax'], policy='one')
    assert (len(st_out) == 2)
    st_out = _post_stmts_preassembly(
        [st1, st2, st3], route, source_apis=['bel', 'biopax'], policy='all')
    assert (len(st_out) == 1)
    st_out = _post_stmts_preassembly(
        [st1, st2, st3], route, source_apis=['bel', 'biopax'], policy='none')
    assert (len(st_out) == 1)


def test_map_grounding():
    route = 'preassembly/map_grounding'
    a = Agent('MEK', db_refs={'TEXT': 'MEK'})
    b = Agent('X', db_refs={'TEXT': 'ERK'})
    st = Activation(a, b)
    st_out = _post_stmts_preassembly([st], route, do_rename=False)
    assert len(st_out) == 1
    assert st_out[0].subj.db_refs.get('FPLX')
    assert st_out[0].obj.db_refs.get('FPLX')
    assert st_out[0].obj.name == 'X'
    st_out = _post_stmts_preassembly([st], route, do_rename=True)
    assert len(st_out) == 1
    assert st_out[0].subj.db_refs.get('FPLX')
    assert st_out[0].obj.db_refs.get('FPLX')
    assert st_out[0].obj.name == 'ERK'


def test_map_grounding_user_map():
    route = 'preassembly/map_grounding'
    gm = {'MEK': {'XXX': 'YYY'}, 'ERK': {'FPLX': 'ERK'}}
    a = Agent('MEK', db_refs={'TEXT': 'MEK'})
    b = Agent('X', db_refs={'TEXT': 'ERK'})
    st = Activation(a, b)
    st_out = _post_stmts_preassembly([st], route, grounding_map=gm,
                                     do_rename=True)
    assert len(st_out) == 1
    assert st_out[0].subj.db_refs.get('XXX') == 'YYY'
    assert st_out[0].obj.db_refs.get('FPLX') == 'ERK'
    assert st_out[0].obj.name == 'ERK'
    gm = {'ERK': {'FPLX': 'ERK_TEST'}}
    st_out = _post_stmts_preassembly([st], route, grounding_map=gm,
                                     grounding_map_policy='extend')
    assert len(st_out) == 1
    assert st_out[0].subj.db_refs.get('FPLX') == 'MEK'
    assert st_out[0].obj.db_refs.get('FPLX') == 'ERK_TEST'
    st_out = _post_stmts_preassembly([st], route, do_rename=True)
    # Make sure the extension to the default grounding map doesn't persist
    assert len(st_out) == 1
    assert st_out[0].subj.db_refs.get('FPLX') == 'MEK'
    assert st_out[0].obj.db_refs.get('FPLX') == 'ERK'
    assert st_out[0].obj.name == 'ERK'


def test_map_sequence():
    route = 'preassembly/map_sequence'
    a = Agent('MAPK1', db_refs={'UP': 'P28482', 'HGNC': '6871'})
    st1 = Phosphorylation(None, a, 'T', '182')
    st2 = Phosphorylation(None, a, 'T', '185')
    st3 = Phosphorylation(None, a, 'Y', '999')
    st_out = _post_stmts_preassembly([st1], route)
    assert len(st_out) == 1, st_out
    assert st_out[0].position == '185'
    st_out = _post_stmts_preassembly([st2], route)
    assert len(st_out) == 1, st_out
    assert st_out[0].position == '185'
    st_out = _post_stmts_preassembly([st3], route)
    assert len(st_out) == 0, st_out


def test_map_sequence_blank_entries():
    """Make sure sites curated as erroneous with no mappings don't
    get treated as valid mappings."""
    route = 'preassembly/map_sequence'
    mapk1 = Agent('MAPK1', db_refs={'UP': 'P28482'})
    rps6 = Agent('RPS6', db_refs={'UP': 'P62753'})
    phos_rps6 = Agent(
        'RPS6', mods=[ModCondition('phosphorylation', 'T', '389')],
        db_refs={'UP': 'P62753'})
    st1 = Phosphorylation(mapk1, rps6, 'T', '389')
    st2 = Phosphorylation(phos_rps6, mapk1, 'T', '185')
    mapped = _post_stmts_preassembly([st1, st2], route)
    assert len(mapped) == 0


def test_filter_by_type():
    route = 'preassembly/filter_by_type'
    st_out = _post_stmts_preassembly(
        [st1, st14], route, stmt_type='Phosphorylation')
    assert len(st_out) == 1


def test_filter_top_level():
    route = 'preassembly/filter_top_level'
    st_out = _post_stmts_preassembly([st14, st15], route)
    assert len(st_out) == 1


def test_filter_no_hypothesis():
    route = 'preassembly/filter_no_hypothesis'
    a = Agent('MAPK1')
    ev1 = Evidence(epistemics={'hypothesis': True})
    ev2 = Evidence(epistemics={'hypothesis': False})
    st1 = Phosphorylation(None, a, evidence=[ev1, ev2])
    st2 = Phosphorylation(None, a, evidence=[ev1, ev1])
    st_out = _post_stmts_preassembly([st1, st2], route)
    assert len(st_out) == 1


def test_filter_no_negated():
    route = 'preassembly/filter_no_negated'
    a = Agent('MAPK1')
    ev1 = Evidence(epistemics={'negated': True})
    ev2 = Evidence(epistemics={'negated': False})
    st1 = Phosphorylation(None, a, evidence=[ev1, ev2])
    st2 = Phosphorylation(None, a, evidence=[ev1, ev1])
    st_out = _post_stmts_preassembly([st1, st2], route)
    assert len(st_out) == 1


def test_belief_cut_plus_filter_top():
    route_belief = 'preassembly/filter_belief'
    route_top_level = 'preassembly/filter_top_level'
    st1 = Phosphorylation(None, Agent('a'))
    st2 = Phosphorylation(Agent('b'), Agent('a'))
    st1.supports = [st2]
    st2.supported_by = [st1]
    st1.belief = 0.9
    st2.belief = 0.1
    st_high_belief = _post_stmts_preassembly(
        [st1, st2], route_belief, belief_cutoff=0.5)
    st_top_level = _post_stmts_preassembly(st_high_belief, route_top_level)
    assert len(st_top_level) == 1


def test_filter_inconsequential_mods():
    route = 'preassembly/filter_inconsequential_mods'
    mc = ModCondition('phosphorylation', None, None, True)
    st1 = Phosphorylation(None, Agent('a'))
    st2 = Phosphorylation(Agent('a', mods=[mc]), Agent('b'))
    st_out = _post_stmts_preassembly([st1, st2], route)
    assert len(st_out) == 1
    whitelist = {'b': [('phosphorylation', None, None)]}
    st_out = _post_stmts_preassembly([st1, st2], route, whitelist=whitelist)
    assert len(st_out) == 2


def test_filter_inconsequential_mods2():
    route = 'preassembly/filter_inconsequential_mods'
    st1 = Phosphorylation(Agent('a'), Agent('b'), 'S', '315')
    whitelist = {'b': [('phosphorylation', 'S', '315')]}
    st_out = _post_stmts_preassembly([st1, st2], route, whitelist=whitelist)
    assert len(st_out) == 1


def test_filter_inconsequential_activities():
    route = 'preassembly/filter_inconsequential_acts'
    st1 = Activation(Agent('a', activity=ActivityCondition('kinase', True)),
                     Agent('b'), 'activity')
    st2 = Activation(Agent('c'), Agent('a'), 'kinase')
    st_out = _post_stmts_preassembly([st1, st2], route)
    assert len(st_out) == 1
    st_out = _post_stmts_preassembly(st_out, route)
    assert len(st_out) == 0


def test_filter_mutation_status():
    route = 'preassembly/filter_mutation_status'
    braf_mut = Agent('BRAF', mutations=MutCondition('600', 'V', 'E'))
    braf_other_mut = Agent('BRAF', mutations=MutCondition('555', 'K', 'G'))
    st1 = Phosphorylation(braf_mut, Agent('a'))
    st2 = Phosphorylation(braf_other_mut, Agent('a'))
    mutations = {'BRAF': [('V', '600', 'E')]}
    deletions = []
    st_out = _post_stmts_preassembly(
        [st1, st2], route, mutations=mutations, deletions=deletions)
    assert len(st_out) == 1
    mutations = {}
    deletions = ['a']
    st_out = _post_stmts_preassembly(
        [st1, st2], route, mutations=mutations, deletions=deletions)
    assert len(st_out) == 0

    # Can we filter statements out based on bound conditions?
    mutations = {'BRAF': [('V', '600', 'E')]}
    deletions = []
    braf_good_bound = deepcopy(braf_mut)
    braf_good_bound.bound_conditions = [BoundCondition(braf_mut)]
    #
    braf_bad_bound = deepcopy(braf_mut)
    braf_bad_bound.bound_conditions = [BoundCondition(braf_other_mut)]
    #
    st3 = Phosphorylation(braf_good_bound, Agent('a'))
    st4 = Phosphorylation(braf_bad_bound, Agent('a'))
    #
    st_out = _post_stmts_preassembly(
        [st3], route, mutations=mutations, deletions=deletions)
    assert len(st_out) == 1
    #
    st_out = _post_stmts_preassembly(
        [st4], route, mutations=mutations, deletions=deletions)
    assert len(st_out) == 0

    # Can we remove bound conditions based on our filter?
    st_out = _post_stmts_preassembly(
        [st3], route, mutations=mutations, deletions=deletions,
        remove_bound=True)
    assert len(st_out[0].enz.bound_conditions) == 1
    #
    st_out = _post_stmts_preassembly(
        [st4], route, mutations=mutations, deletions=deletions,
        remove_bound=True)
    assert len(st_out[0].enz.bound_conditions) == 0


def test_rename_db_ref():
    route = 'preassembly/rename_db_ref'
    x = Agent('X', db_refs={'BE': 'X'})
    y = Agent('Y', db_refs={'FPLX': 'Y'})
    st1 = Phosphorylation(x, y)
    st_out = _post_stmts_preassembly([st1], route, ns_from='BE', ns_to='FPLX')
    assert len(st_out) == 1
    assert st_out[0].enz.db_refs.get('FPLX') == 'X'
    assert 'BE' not in st_out[0].enz.db_refs
    assert st_out[0].sub.db_refs.get('FPLX') == 'Y'


def test_filter_concept_names():
    route = 'preassembly/filter_concept_names'
    stmts = [
        Influence(Event(Concept('a')), Event(Concept('b'))),
        Influence(Event(Concept('a')), Event(Concept('c'))),
        Influence(Event(Concept('a')), Event(Concept('d'))),
        Influence(Event(Concept('c')), Event(Concept('d')))
        ]
    st_out = _post_stmts_preassembly(
        stmts, route, name_list=['a'], policy='one')
    assert len(st_out) == 3, st_out
    st_out = _post_stmts_preassembly(
        stmts, route, name_list=['a', 'b', 'c'], policy='all')
    assert len(st_out) == 2, st_out
    st_out = _post_stmts_preassembly(
        stmts, route, name_list=['a', 'd'], policy='one')
    assert len(st_out) == 4, st_out
    st_out = _post_stmts_preassembly(
        stmts, route, name_list=['a', 'b'], policy='one', invert=True)
    assert len(st_out) == 1, st_out


def test_filter_namespace_concepts_simple():
    def make_statement(a, b):
        return Influence(Event(Concept(a, db_refs={'TEXT': a})),
                         Event(Concept(b, db_refs={'TEXT': b})))
    route = 'preassembly/filter_by_db_refs'
    stmts = [make_statement('education', 'thinking'),
             make_statement('doubt', 'government')]
    fs = _post_stmts_preassembly(
        stmts, route, namespace='TEXT', values=['education'], policy='one')
    assert len(fs) == 1
    assert stmts[0].equals(fs[0])
    fs = _post_stmts_preassembly(
        stmts, route, namespace='TEXT', values=['education'], policy='one',
        invert=True)
    assert len(fs) == 2
    assert stmts[0].equals(fs[0])
    assert stmts[1].equals(fs[1])
    fs = _post_stmts_preassembly(
        stmts, route, namespace='TEXT', values=['education'], policy='all',
        invert=True)
    assert len(fs) == 1
    assert stmts[1].equals(fs[0])
    fs = _post_stmts_preassembly(
        stmts, route, namespace='TEXT', values=['education'], policy='all')
    assert not fs, fs


def test_filter_namespace_concepts_list():
    def make_statement(a, b):
        return Influence(Event(Concept(a, db_refs={'UN': [(a, 1.0)]})),
                         Event(Concept(b, db_refs={'UN': [(b, 1.0)]})))
    route = 'preassembly/filter_by_db_refs'
    stmts = [make_statement('UN/entities/human/education',
                            'UN/entities/human/food/food_security'),
             make_statement('UN/entities/human/fishery',
                            'UN/entities/human/government')]
    fs = _post_stmts_preassembly(
        stmts, route, namespace='UN', values=['education'], policy='one',
        match_suffix=True)
    assert len(fs) == 1
    assert stmts[0].equals(fs[0])
    fs = _post_stmts_preassembly(
        stmts, route, namespace='UN', values=['education', 'fishery'],
        policy='one',  match_suffix=True)
    assert len(fs) == 2
    assert stmts[0].equals(fs[0])
    assert stmts[1].equals(fs[1])
    fs = _post_stmts_preassembly(
        stmts, route, namespace='UN', values=['fishery', 'government'],
        policy='all',  match_suffix=True)
    assert len(fs) == 1
    assert stmts[1].equals(fs[0])


def test_merge_groundings():
    route_preassembly = 'preassembly/run_preassembly'
    route_groundings = 'preassembly/merge_groundings'
    refs1 = {'UN': [('x', 0.8), ('y', 0.7)],
             'B': 'x',
             'C': 'y'}
    refs2 = {'UN': [('x', 0.9), ('y', 0.6), ('z', 0.5)],
             'B': 'x',
             'D': 'z'}
    stmts = [Influence(Event(Concept('a', db_refs=refs1)),
                       Event(Concept('b', db_refs=refs2)),
                       evidence=[Evidence(source_api='eidos', text='1')]),
             Influence(Event(Concept('a', db_refs=refs2)),
                       Event(Concept('b', db_refs=refs1)),
                       evidence=[Evidence(source_api='eidos', text='2')])]
    stmts = _post_stmts_preassembly(stmts, route_preassembly)
    assert len(stmts) == 1
    stmts = _post_stmts_preassembly(stmts, route_groundings)
    assert stmts[0].subj.concept.db_refs == \
        {'UN': [('x', 0.9), ('y', 0.7), ('z', 0.5)],
         'B': 'x', 'C': 'y', 'D': 'z'}, \
        stmts[0].subj.db_refs
    assert stmts[0].obj.concept.db_refs == stmts[0].subj.concept.db_refs


def test_merge_deltas():
    def add_annots(stmt):
        for ev in stmt.evidence:
            ev.annotations['subj_adjectives'] = stmt.subj.delta.adjectives
            ev.annotations['obj_adjectives'] = stmt.obj.delta.adjectives
            ev.annotations['subj_polarity'] = stmt.subj.delta.polarity
            ev.annotations['obj_polarity'] = stmt.obj.delta.polarity
        return stmt

    d1 = QualitativeDelta(polarity=1, adjectives=['a', 'b', 'c'])
    d2 = QualitativeDelta(polarity=-1, adjectives=None)
    d3 = QualitativeDelta(polarity=1, adjectives=['g'])
    d4 = QualitativeDelta(polarity=-1, adjectives=['d', 'e', 'f'])
    d5 = QualitativeDelta(polarity=None, adjectives=['d'])
    d6 = QualitativeDelta(polarity=None, adjectives=None)
    d7 = QualitativeDelta(polarity=1, adjectives=None)

    def make_ev(name, delta):
        return Event(Concept(name), delta=delta)

    route_preassembly = 'preassembly/run_preassembly'
    route_deltas = 'preassembly/merge_deltas'

    stmts = [add_annots(Influence(make_ev('a', sd), make_ev('b', od),
                                  evidence=[Evidence(source_api='eidos',
                                                     text='%d' % idx)]))
             for idx, (sd, od) in enumerate([(d1, d2), (d3, d4)])]
    stmts = _post_stmts_preassembly(
        stmts, route_preassembly, return_toplevel=True)
    stmts = _post_stmts_preassembly(stmts, route_deltas)
    assert stmts[0].subj.delta.polarity == 1, stmts[0].subj.delta
    assert stmts[0].obj.delta.polarity == -1, stmts[0].obj.delta
    assert set(stmts[0].subj.delta.adjectives) == {'a', 'b', 'c', 'g'}, \
        stmts[0].subj.delta
    assert set(stmts[0].obj.delta.adjectives) == {'d', 'e', 'f'}, \
        stmts[0].obj.delta

    stmts = [add_annots(Influence(make_ev('a', sd), make_ev('b', od),
                                  evidence=[Evidence(source_api='eidos',
                                                     text='%d' % idx)]))
             for idx, (sd, od) in enumerate([(d1, d5), (d6, d7), (d6, d7)])]
    stmts = _post_stmts_preassembly(
        stmts, route_preassembly, return_toplevel=True)
    stmts = _post_stmts_preassembly(stmts, route_deltas)
    assert stmts[0].subj.delta.polarity is None, stmts[0].subj.delta
    assert stmts[0].obj.delta.polarity == 1, stmts[0].obj.delta
    assert set(stmts[0].subj.delta.adjectives) == {'a', 'b', 'c'}, \
        stmts[0].subj.delta
    assert set(stmts[0].obj.delta.adjectives) == {'d'}, \
        stmts[0].obj.delta


def test_preassemble_flatten():
    route = 'preassembly/run_preassembly'
    st_out = _post_stmts_preassembly(
        [st1, st3, st5, st6], route, flatten_evidence=False)
    assert len(st_out[0].evidence) == 1
    assert len(st_out[1].evidence) == 1
    st_out = _post_stmts_preassembly(
        [st1, st3, st5, st6], route, flatten_evidence=True,
        flatten_evidence_collect_from='supported_by')
    assert len(st_out[0].evidence) == 2
    assert len(st_out[1].evidence) == 2
    st_out = _post_stmts_preassembly(
        [st1, st3, st5, st6], route, flatten_evidence=True,
        flatten_evidence_collect_from='supports')
    assert len(st_out[0].evidence) == 1
    assert len(st_out[1].evidence) == 1


def test_filter_by_curation():
    route = 'preassembly/filter_by_curation'
    new_st1 = deepcopy(st1)
    new_ev = Evidence(text='a -> b', source_api='new')
    new_st1.evidence.append(new_ev)
    stmts_in = [new_st1, st2, st3]
    assert len(new_st1.evidence) == 2
    assert all(st.belief != 1 for st in stmts_in)
    cur1 = {'pa_hash': new_st1.get_hash(),
            'source_hash': new_st1.evidence[0].get_source_hash(),
            'tag': 'grounding'}
    cur2 = {'pa_hash': new_st1.get_hash(),
            'source_hash': new_st1.evidence[1].get_source_hash(),
            'tag': 'wrong_relation'}
    cur3 = {'pa_hash': new_st1.get_hash(),
            'source_hash': new_st1.evidence[0].get_source_hash(),
            'tag': 'correct'}
    cur4 = {'pa_hash': st2.get_hash(),
            'source_hash': st2.evidence[0].get_source_hash(),
            'tag': 'correct'}
    # With 'any' policy it is enough to have one incorrect curation
    any_incorrect_one_cur = _post_stmts_preassembly(
        stmts_in, route, curations=[cur1], incorrect_policy='any')
    assert len(any_incorrect_one_cur) == 2
    assert not any(
        [new_st1.get_hash() == st.get_hash() for st in any_incorrect_one_cur])
    # With 'all' policy all evidences have to be curated
    all_incorrect_one_cur = _post_stmts_preassembly(
        stmts_in, route, curations=[cur1], incorrect_policy='all')
    assert len(all_incorrect_one_cur) == 3, len(all_incorrect_one_cur)
    assert any(
        [new_st1.get_hash() == st.get_hash() for st in all_incorrect_one_cur])
    all_incorrect_two_cur = _post_stmts_preassembly(
        stmts_in, route, curations=[cur1, cur2], incorrect_policy='all')
    assert len(all_incorrect_two_cur) == 2
    assert not any(
        [new_st1.get_hash() == st.get_hash() for st in all_incorrect_two_cur])
    # Correct curation cancels out incorrect
    correct_incorrect = _post_stmts_preassembly(
        stmts_in, route, curations=[cur1, cur2, cur3, cur4],
        incorrect_policy='all', update_belief=False)
    assert len(correct_incorrect) == 3, len(correct_incorrect)
    assert any(
        [new_st1.get_hash() == st.get_hash() for st in correct_incorrect])
    assert all(st.belief != 1 for st in correct_incorrect)
    # Optionally update belief to 1 for correct curation
    new_belief = _post_stmts_preassembly(
        stmts_in, route, curations=[cur1, cur2, cur3, cur4],
        incorrect_policy='all', update_belief=True)
    assert new_belief[0].belief == 1
    assert new_belief[1].belief == 1
    assert new_belief[2].belief == 0.7


def test_eidos_ungrounded():
    route = 'preassembly/filter_grounded_only'
    a = Agent('x', db_refs={'TEXT': 'x', 'TEXT_NORM': 'y'})
    b = Agent('x', db_refs={'TEXT': 'x', })
    c = Agent('x', db_refs={'TEXT': 'x', 'GO': 'GO:1234'})
    stmts = [Activation(a, b),
             Activation(a, c),
             Activation(b, c),
             Activation(c, c)]
    st_out = _post_stmts_preassembly(stmts, route)
    assert len(st_out) == 1


def test_responsive():
    res = _call_api('get', '/')
    assert res.data.startswith(
        b'<!DOCTYPE html>\n<html>\n<head>\n    <title>INDRA REST API</title>'), \
        "Unexpected content: %s" % res.data


def test_options():
    res = _call_api('options', 'reach/process_text')
    res_json = json.loads(res.get_data())
    assert res_json == {}, \
        "Unexpected content: %s" % res_json


@attr('notravis')
def test_trips_process_text():
    res = _call_api('post', 'trips/process_text',
                    json={'text': 'MEK phosphorylates ERK.'})
    res_json = json.loads(res.get_data())
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert len(stmts) == 1, len(stmts)
    stmt = stmts[0]
    assert isinstance(stmt, Phosphorylation), type(stmt)
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.sub.name == 'ERK', stmt.sub


def test_trips_process_xml():
    test_ekb_path = path.join(HERE, 'trips_ekbs',
                              'MEK_increases_the_phosphorylation_of_ERK.ekb')
    with open(test_ekb_path, 'r') as f:
        xml_str = f.read()
    res = _call_api('post', 'trips/process_xml', json={'xml_str': xml_str})
    res_json = json.loads(res.get_data())
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert len(stmts) == 1, len(stmts)
    stmt = stmts[0]
    assert isinstance(stmt, Phosphorylation), type(stmt)
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.sub.name == 'ERK', stmt.sub


@attr('notravis')
def test_reach_process_text():
    res = _call_api('post', 'reach/process_text',
                    json={'text': 'MEK phosphorylates ERK.'})
    res_json = json.loads(res.get_data())
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert len(stmts) == 1, len(stmts)
    stmt = stmts[0]
    assert isinstance(stmt, Phosphorylation), type(stmt)
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.sub.name == 'ERK', stmt.sub


def test_reach_process_json():
    test_file = path.join(HERE, 'reach_act_amt.json')
    with open(test_file, 'rb') as fh:
        json_str = fh.read().decode('utf-8')
    res = _call_api('post', 'reach/process_json',
                    json={'json': json_str})
    res_json = json.loads(res.get_data())
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert len(stmts) == 1
    assert isinstance(stmts[0], IncreaseAmount)
    assert stmts[0].subj is not None
    assert stmts[0].obj is not None


def test_reach_process_pmc():
    res = _call_api('post', 'reach/process_pmc', json={'pmcid': 'PMC4338247'})
    res_json = json.loads(res.get_data())
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert stmts
    assert stmts[0].evidence[0].pmid is not None


def test_cwms_process_text():
    res = _call_api('post', '/cwms/process_text',
                    json={'text': 'Hunger causes displacement.'})
    res_json = json.loads(res.get_data())
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert len(stmts) == 1


def test_hume_process_jsonld():
    from indra.tests.test_hume import test_file_new_simple
    with open(test_file_new_simple, 'r') as fh:
        test_jsonld = fh.read()
    res = _call_api('post', '/hume/process_jsonld',
                    json={'jsonld': test_jsonld})
    res_json = json.loads(res.get_data())
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert len(stmts) == 1


def test_eidos_json():
    from indra.tests.test_eidos import test_jsonld
    with open(test_jsonld, 'r') as fh:
        test_json = fh.read()
    res = _call_api('post', '/eidos/process_jsonld',
                    json={'jsonld': test_json})
    res_json = json.loads(res.get_data())
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert len(stmts) == 1


STMT_JSON = {
    "id": "acc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
    "type": "Complex",
    "members": [
        {"db_refs": {"TEXT": "MEK", "FPLX": "MEK"}, "name": "MEK"},
        {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"}
    ],
    "sbo": "http://identifiers.org/sbo/SBO:0000526",
    "evidence": [{"text": "MEK binds ERK", "source_api": "trips"}]
}


def test_assemblers_cyjs():
    res = _call_api('post', 'assemblers/cyjs',
                    json={'statements': [STMT_JSON]})
    res_json = json.loads(res.get_data())
    print(res_json)
    assert len(res_json['edges']) == 1, len(res_json['edges'])
    assert len(res_json['nodes']) == 2, len(res_json['nodes'])
    return


def test_assemblers_pysb_no_format():
    res = _call_api('post', 'assemblers/pysb',
                    json={'statements': [STMT_JSON]})
    res_json = json.loads(res.get_data())
    assert 'model' in res_json.keys()
    assert res_json['model'] is not None
    assert 'MEK' in res_json['model'], res_json['model']
    return


def test_assemblers_pysb_kappa_img_format():
    for exp_format in ['kappa_im', 'kappa_cm']:
        print("Testing", exp_format)
        res = _call_api('post', 'assemblers/pysb', json={
            'statements': [STMT_JSON], 'export_format': exp_format})
        res_json = json.loads(res.get_data())
        assert 'image' in res_json.keys()
        assert 'image' in res_json.keys()
        assert res_json['image'] is not None
    return


def test_assemblers_pysb_kappa_other_formats():
    # All the formats defined in PysbAssembler.export_model doc string.
    formats = ['bngl', 'kappa', 'sbml', 'matlab', 'mathematica',
               'potterswheel']
    for exp_format in formats:
        print("Testing", exp_format)
        res = _call_api('post', 'assemblers/pysb', json={
            'statements': [STMT_JSON], 'export_format': exp_format})
        res_json = json.loads(res.get_data())
        assert 'model' in res_json.keys()
        assert res_json['model'] is not None
        assert 'MEK' in res_json['model'], res_json['model']
    return


def test_assemblers_cx():
    res = _call_api('post', 'assemblers/cx', json={'statements': [STMT_JSON]})
    res_json = json.loads(res.get_data())
    assert 'model' in res_json.keys()
    assert res_json['model'] is not None
    assert 'MEK' in res_json['model'], res_json['model']


def test_assemblers_graph():
    res = _call_api('post', 'assemblers/graph',
                    json={'statements': [STMT_JSON]})
    res_json = json.loads(res.get_data())
    assert 'model' in res_json.keys()
    assert res_json['model'] is not None
    assert 'MEK' in res_json['model'], res_json['model']


def test_assemblers_english():
    res = _call_api('post', 'assemblers/english',
                    json={'statements': [STMT_JSON]})
    res_json = json.loads(res.get_data())
    assert 'sentences' in res_json.keys()
    assert len(res_json['sentences']) == 1, len(res_json['sentences'])
    sentence = list(res_json['sentences'].values())[0]
    assert 'MEK' in sentence, sentence


def test_assemblers_loopy():
    stmt_jsons = [{
            "id": "acc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
            "type": "Phosphorylation",
            "enz": {"db_refs": {"TEXT": "MEK", "FPLX": "MEK"}, "name": "MEK"},
            "sub": {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"},
            "sbo": "http://identifiers.org/sbo/SBO:0000526",
            "evidence": [
                {"text": "MEK phosphorylates ERK", "source_api": "trips"}]
        },
        {
            "id": "bcc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
            "type": "Activation",
            "subj": {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"},
            "obj": {"db_refs": {"TEXT": "EGFR", "HGNC": "3236"},
                    "name": "EGFR"},
            "sbo": "http://identifiers.org/sbo/SBO:0000526",
            "evidence": [{"text": "ERK activates EGFR", "source_api": "trips"}]
        }
    ]
    res = _call_api('post', 'assemblers/sif/loopy',
                    json={'statements': stmt_jsons})
    res_json = json.loads(res.get_data())
    assert 'loopy_url' in res_json.keys()
    assert "ERK" in res_json['loopy_url']


def test_pipeline():
    p = [{'function': 'filter_grounded_only'},
         {'function': 'run_preassembly',
          'kwargs': {'return_toplevel': False}}]
    res = _call_api('post', 'preassembly/pipeline',
                    json={'statements': [STMT_JSON], 'pipeline': p})
    res_json = json.loads(res.get_data())
    assert 'statements' in res_json
    assert len(res_json['statements']) == 1
