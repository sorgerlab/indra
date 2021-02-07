from indra.statements import *
from indra.tools.adaptive_assembly import AdaptiveAssembler
from indra.preassembler import OntologyRefinementFilter, \
    RefinementConfirmationFilter
from indra.ontology.bio import bio_ontology


def test_adaptive_assembly():
    erk = Agent('ERK', db_refs={'FPLX': 'ERK'})
    mek = Agent('MEK', db_refs={'FPLX': 'MEK'})
    mapk = Agent('MAPK', db_refs={'FPLX': 'MAPK'})

    stmts = [
        Phosphorylation(mek, erk),
        Phosphorylation(mek, mapk, 'T'),
        Phosphorylation(mek, erk, 'T')
    ]
    hashes = [stmt.get_hash() for stmt in stmts]
    filters = [
        OntologyRefinementFilter(ontology=bio_ontology),
        RefinementConfirmationFilter(ontology=bio_ontology),
    ]
    aa = AdaptiveAssembler(stmts, filters=filters)

    all_refinements = aa.get_all_refinements()
    assert set(all_refinements) == {(hashes[2], hashes[0]),
                                    (hashes[2], hashes[1])}

    test_stmt = Phosphorylation(mek, mapk)
    test_refinements = aa.get_refinements(test_stmt)
    assert False, test_refinements