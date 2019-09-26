import unittest
from indra.statements import Influence, Concept, Event
from indra.preassembler.ontology_mapper import OntologyMapper, wm_ontomap


def test_map():
    c1 = Concept('x', db_refs={'UN': [('entities/x', 1.0)]})
    c2 = Concept('y', db_refs={'HUME': [('entities/y', 1.0)]})
    c3 = Concept('z')
    stmts = [Influence(Event(c1), Event(c3)),
             Influence(Event(c2), Event(c3))]
    mappings = [(('UN', 'entities/x'), ('HUME', 'entities/y'))]
    om = OntologyMapper(stmts, mappings=mappings)
    om.map_statements()
    assert len(om.statements) == 2
    assert om.statements[0].subj.concept.db_refs['HUME'] == \
        [('entities/y', 1.0)], om.statements[0].subj.concept.db_refs
    assert om.statements[1].subj.concept.db_refs['UN'] == \
        [('entities/x', 1.0)], om.statements[1].subj.concept.db_refs


def test_wm_map():
    famine = 'wm/concept/causal_factor/condition/famine'
    edu = 'wm/concept/causal_factor/social_and_political/education/education'
    c1 = Concept('x', db_refs={'WM': [(famine, 1.0)]})
    c2 = Concept('y', db_refs={'WM': [(edu, 1.0)]})
    stmts = [Influence(Event(c1), Event(c2))]
    om = OntologyMapper(stmts, wm_ontomap, symmetric=False)
    om.map_statements()
    stmt = om.statements[0]
    assert 'SOFIA' in stmt.subj.concept.db_refs
    assert 'SOFIA' in stmt.obj.concept.db_refs

    # Test the previously problematic famine case
    c3 = Concept('z', db_refs={'SOFIA': 'events1/Health/Famine'})

    # Unscored mapping
    om = OntologyMapper([Event(c3)], wm_ontomap, symmetric=False, scored=False)
    om.map_statements()
    stmt = om.statements[0]
    assert stmt.concept.db_refs['WM'] == [(famine, 1.0)], \
        stmt.concept.db_refs['WM']

    c3 = Concept('z', db_refs={'SOFIA': 'events1/Health/Famine'})
    om = OntologyMapper([Event(c3)], wm_ontomap, symmetric=False, scored=True)
    om.map_statements()
    stmt = om.statements[0]
    assert stmt.concept.db_refs['WM'] == [(famine, 0.79034126)], \
        stmt.concept.db_refs['WM']
