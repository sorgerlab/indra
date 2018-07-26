import unittest
from indra.statements import Influence, Concept
from indra.preassembler.ontology_mapper import OntologyMapper, wm_ontomap


def test_map():
    c1 = Concept('x', db_refs={'UN': [('entities/x', 1.0)]})
    c2 = Concept('y', db_refs={'BBN': [('entities/y', 1.0)]})
    c3 = Concept('z')
    stmts = [Influence(c1, c3), Influence(c2, c3)]
    om = OntologyMapper(stmts)
    om.map_statements()
    assert len(om.statements) == 2
    assert om.statements[0].subj.db_refs['BBN'] == [('entities/y', 1.0)] \
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
    assert 'BBN' in stmt.subj.db_refs
    assert 'BBN' in stmt.obj.db_refs
    assert 'SOFIA' in stmt.subj.db_refs
    assert 'SOFIA' in stmt.obj.db_refs
