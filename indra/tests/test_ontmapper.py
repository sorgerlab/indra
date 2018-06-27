from indra.statements import Influence, Concept
from indra.preassembler.ontology_mapper import OntologyMapper, wm_ontomap


def test_map():
    c1 = Concept('x', db_refs={'UN': [('entities/x', 1.0)]})
    c2 = Concept('y', db_refs={'BBN': 'entities/y'})
    c3 = Concept('z')
    stmts = [Influence(c1, c3), Influence(c2, c3)]
    om = OntologyMapper(stmts)
    om.map_statements()
    assert len(om.statements) == 2
    assert om.statements[0].subj.db_refs['BBN'] == 'entities/y', \
        om.statements[0].subj.db_refs
    assert om.statements[1].subj.db_refs['UN'] == [('entities/x', 1.0)], \
        om.statements[1].subj.db_refs


def test_wm_map():
    c1 = Concept('x', db_refs={'UN': [('')]})
    om = OntologyMapper(stmts, wm_ontomap)