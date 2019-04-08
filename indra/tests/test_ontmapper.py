import unittest
from indra.statements import Influence, Concept, Event
from indra.preassembler.ontology_mapper import OntologyMapper, wm_ontomap


def test_map():
    c1 = Concept('x', db_refs={'UN': [('entities/x', 1.0)]})
    c2 = Concept('y', db_refs={'HUME': [('entities/y', 1.0)]})
    c3 = Concept('z')
    stmts = [Influence(Event(c1), Event(c3)),
             Influence(Event(c2), Event(c3))]
    om = OntologyMapper(stmts)
    om.map_statements()
    assert len(om.statements) == 2
    assert om.statements[0].subj.concept.db_refs['HUME'] == \
        [('entities/y', 1.0)], om.statements[0].subj.concept.db_refs
    assert om.statements[1].subj.concept.db_refs['UN'] == \
        [('entities/x', 1.0)], om.statements[1].subj.concept.db_refs


def test_wm_map():
    c1 = Concept('x', db_refs={'UN': [('UN/events/human/famine', 1.0)]})
    c2 = Concept('y', db_refs={'UN': [('UN/entities/human/education', 1.0)]})
    stmts = [Influence(Event(c1), Event(c2))]
    om = OntologyMapper(stmts, wm_ontomap, symmetric=False)
    om.map_statements()
    stmt = om.statements[0]
    assert 'HUME' in stmt.subj.concept.db_refs
    assert 'HUME' in stmt.obj.concept.db_refs
    assert 'SOFIA' in stmt.subj.concept.db_refs
    assert 'SOFIA' in stmt.obj.concept.db_refs

    # Test the previously problematic famine case
    c3 = Concept('z', db_refs={'SOFIA': 'Health/Famine'})
    c4 = Concept('a', db_refs={'HUME': [('event/healthcare/famine', 1.0)]})
    stmts = [Influence(Event(c4), Event(c3))]

    # Unscored mapping
    om = OntologyMapper(stmts, wm_ontomap, symmetric=False, scored=False)
    om.map_statements()
    stmt = om.statements[0]
    assert stmt.obj.concept.db_refs['UN'] == [('UN/events/human/famine',
                                               1.0)], \
        stmt.obj.concept.db_refs['UN']
    assert stmt.subj.concept.db_refs['UN'] == [('UN/events/human/famine',
                                                1.0)], \
        stmt.subj.concept.db_refs['UN']

    # Scored mapping
    c3 = Concept('z', db_refs={'SOFIA': 'Health/Famine'})
    c4 = Concept('a', db_refs={'HUME': [('event/healthcare/famine', 1.0)]})
    stmts = [Influence(Event(c4), Event(c3))]
    om = OntologyMapper(stmts, wm_ontomap, symmetric=False, scored=True)
    om.map_statements()
    stmt = om.statements[0]
    assert stmt.obj.concept.db_refs['UN'] == [('UN/events/human/famine',
                                               0.81851065)], \
        stmt.obj.concept.db_refs['UN']
    assert stmt.subj.concept.db_refs['UN'] == [('UN/events/human/famine',
                                                1.0)], \
        stmt.subj.concept.db_refs['UN']
