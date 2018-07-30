import unittest
from indra.statements import Influence, Concept
from indra.preassembler.ontology_mapper import OntologyMapper, wm_ontomap


def test_map():
    c1 = Concept('x', db_refs={'UN': [('entities/x', 1.0)]})
    c2 = Concept('y', db_refs={'HUME': [('entities/y', 1.0)]})
    c3 = Concept('z')
    stmts = [Influence(c1, c3), Influence(c2, c3)]
    om = OntologyMapper(stmts)
    om.map_statements()
    assert len(om.statements) == 2
    assert om.statements[0].subj.db_refs['HUME'] == [('entities/y', 1.0)], \
        om.statements[0].subj.db_refs
    assert om.statements[1].subj.db_refs['UN'] == [('entities/x', 1.0)], \
        om.statements[1].subj.db_refs


@unittest.skip('Mapping file not in repo')
def test_wm_map():
    c1 = Concept('x', db_refs={'UN': [('UN/properties/price', 1.0)]})
    c2 = Concept('y', db_refs={'UN': [('UN/entities/human/education', 1.0)]})
    stmts = [Influence(c1, c2)]
    om = OntologyMapper(stmts, wm_ontomap, symmetric=False)
    om.map_statements()
    stmt = om.statements[0]
    assert 'HUME' in stmt.subj.db_refs
    assert 'HUME' in stmt.obj.db_refs
    assert 'SOFIA' in stmt.subj.db_refs
    assert 'SOFIA' in stmt.obj.db_refs

    # Test the previously problematic famine case
    c3 = Concept('z', db_refs={'SOFIA': 'Health/Famine'})
    c4 = Concept('a', db_refs={'HUME': [('event/healthcare/famine', 1.0)]})
    stmts = [Influence(c4, c3)]

    # Unscored mapping
    om = OntologyMapper(stmts, wm_ontomap, symmetric=False, scored=False)
    om.map_statements()
    stmt = om.statements[0]
    assert stmt.obj.db_refs['UN'] == [('UN/events/human/famine', 1.0)], \
        stmt.obj.db_refs['UN']
    assert stmt.subj.db_refs['UN'] == [('UN/events/human/famine', 1.0)], \
        stmt.subj.db_refs['UN']

    # Scored mapping
    c3 = Concept('z', db_refs={'SOFIA': 'Health/Famine'})
    c4 = Concept('a', db_refs={'HUME': [('event/healthcare/famine', 1.0)]})
    stmts = [Influence(c4, c3)]
    om = OntologyMapper(stmts, wm_ontomap, symmetric=False, scored=True)
    om.map_statements()
    stmt = om.statements[0]
    assert stmt.obj.db_refs['UN'] == [('UN/events/human/famine',
                                       1.036856450549298)], \
        stmt.obj.db_refs['UN']
    assert stmt.subj.db_refs['UN'] == [('UN/events/human/famine',
                                        0.9465582395117739)], \
        stmt.subj.db_refs['UN']
