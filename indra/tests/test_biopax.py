import os
from collections import defaultdict
from indra.sources import biopax
from indra.statements import *
import indra.sources.biopax.processor as bpc
from indra.util import unicode_strs
import pytest

model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'biopax_test.owl')

bp = biopax.process_owl(model_path)

stmts_by_source_id = defaultdict(set)
for stmt in bp.statements:
    for ev in stmt.evidence:
        stmts_by_source_id[ev.source_id.split('/')[-1]].add(stmt)


def test_listify():
    assert bpc._listify(1) == [1]
    assert bpc._listify([1,2] == [1,2])
    assert bpc._listify([1] == [1])


def test_protein_family_agent():
    bpe = bp.model.objects['Protein_da79d1a005a8eb259b0c09278ae9230e']
    agents = bp._get_agents_from_entity(bpe)
    assert len(agents) == 2
    assert {a.name for a in agents} == {'MAPK1', 'MAPK3'}

    mapk3 = agents[0]
    assert mapk3.name == 'MAPK3'
    assert len(mapk3.mods) == 2
    assert mapk3.mods[0].position == '202'
    assert mapk3.db_refs['UP'] == 'P27361', mapk3.db_refs
    assert mapk3.db_refs['HGNC'] == '6877', mapk3.db_refs


def test_phosphorylation_extraction():
    cat = 'Catalysis_12411c77e2fd2252f2a8b52bbee6eeb9'
    stmts = list(stmts_by_source_id[cat])
    assert len(stmts) == 2
    assert all(isinstance(stmt, Phosphorylation) for stmt in stmts)
    assert all(stmt.enz.name == 'BRAF' for stmt in stmts)
    assert all(stmt.sub.name == 'MAP2K1' for stmt in stmts)
    assert {stmt.position for stmt in stmts} == {'218', '222'}


def test_active_form_extraction():
    br = 'BiochemicalReaction_5271ab5d568f463bffa8bb378e8aa257'
    stmts = list(stmts_by_source_id[br])
    assert len(stmts) == 2, stmts
    stmt = stmts[0]
    assert isinstance(stmt, ActiveForm)
    assert stmt.agent.name in {'MAPK1', 'MAPK3'}
    assert len(stmt.agent.mods) == 2


def test_amount_regulation_extraction():
    br = 'TemplateReactionRegulation_77553a64e29e82b517d0be230363b757'
    stmts = list(stmts_by_source_id[br])
    assert len(stmts) == 1, stmts
    stmt = stmts[0]
    assert isinstance(stmt, IncreaseAmount)
    assert stmt.subj.name == 'reactive oxygen species'
    assert stmt.subj.db_refs['CHEBI'] == 'CHEBI:26523'
    assert stmt.obj.name == 'TP53'


def test_activity_regulation_extraction():
    br = 'Catalysis_5710598d0317be9d4742b469d1322a48'
    stmts = list(stmts_by_source_id[br])
    assert len(stmts) == 3, stmts
    assert {s.subj.name for s in stmts} == {'HRAS', 'NRAS', 'KRAS'}
    assert isinstance(stmts[0], Activation)
    assert all(s.obj.name == 'BRAF' for s in stmts)


def test_chebi_grounding_extraction():
    bpe = 'SmallMolecule_49d78305d95647ad81961ec7f6189821'
    sm = bp.model.objects[bpe]
    agents = bp._get_agents_from_singular_entity(sm)
    assert len(agents) == 1
    assert agents[0].name == 'GTP'
    assert agents[0].db_refs['CHEBI'] == 'CHEBI:15996'


@pytest.mark.webservice
@pytest.mark.slow
def test_pathsfromto():
    bp = biopax.process_pc_pathsfromto(['MAP2K1'], ['MAPK1'])
    assert_pmids(bp.statements)
    assert_source_sub_id(bp.statements)
    assert unicode_strs(bp.statements)
    num_unique = len({s.get_hash(shallow=False) for s in bp.statements})
    assert len(bp.statements) == num_unique


def assert_pmids(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid is not None:
                assert ev.pmid.isdigit()


def assert_source_sub_id(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            assert 'source_sub_id' in ev.annotations
            assert ev.annotations['source_sub_id']


def test_valid_agent():
    agent = bpc.get_standard_agent('x', {'HGNC': '1097', 'EGID': '---'})
    assert agent.name == 'BRAF'
    assert agent.db_refs.get('EGID') != '---'
    assert agent.db_refs['HGNC'] == '1097'
