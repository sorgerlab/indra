from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra.preassembler import Preassembler, render_stmt_graph, \
                               flatten_evidence, flatten_stmts
from indra import trips
from indra.statements import Agent, Phosphorylation, BoundCondition, \
                             Dephosphorylation, Evidence, ModCondition, \
                             ActiveForm, MutCondition, Complex, Translocation
from indra.preassembler.hierarchy_manager import hierarchies

"""
def test_from_text():
    sentences = ['Src phosphorylates Ras, bound to GTP, at Tyr32.', 
                'Src phosphorylates NRAS at Tyr32.',
                'Src phosphorylates NRAS that is bound to GTP.']
    pa = Preassembler()
    for s in sentences:
        tp = trips.process_text(s)
        pa.add_statements(tp.statements)
"""

def test_duplicates():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FA': '03663'})
    st1 = Phosphorylation(src, ras)
    st2 = Phosphorylation(src, ras)
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    pa.combine_duplicates()
    assert(len(pa.unique_stmts) == 1)

def test_duplicates_copy():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FA': '03663'})
    st1 = Phosphorylation(src, ras, evidence=[Evidence(text='Text 1')])
    st2 = Phosphorylation(src, ras, evidence=[Evidence(text='Text 2')])
    stmts = [st1, st2]
    pa = Preassembler(hierarchies, stmts=stmts)
    pa.combine_duplicates()
    assert(len(pa.unique_stmts) == 1)
    assert(len(stmts) == 2)
    assert(len(stmts[0].evidence) == 1)
    assert(len(stmts[1].evidence) == 1)

def test_duplicates_sorting():
    mc = ModCondition('phosphorylation')
    map2k1_1 = Agent('MAP2K1', mods=[mc])
    mc1 = ModCondition('phosphorylation', 'serine', '218')
    mc2 = ModCondition('phosphorylation', 'serine', '222')
    mc3 = ModCondition('phosphorylation', 'serine', '298')
    map2k1_2 = Agent('MAP2K1', mods=[mc1, mc2, mc3])
    mapk3 = Agent('MAPK3')
    #ras = Agent('MAPK3', db_refs = {'FA': '03663'})
    #nras = Agent('NRAS', db_refs = {'FA': '03663'})
    st1 = Phosphorylation(map2k1_1, mapk3, position='218')
    st2 = Phosphorylation(map2k1_2, mapk3)
    st3 = Phosphorylation(map2k1_1, mapk3, position='218')
    stmts = [st1, st2, st3]
    pa = Preassembler(hierarchies, stmts=stmts)
    pa.combine_duplicates()
    assert(len(pa.unique_stmts) == 2)

def test_combine_duplicates():
    raf = Agent('RAF1')
    mek = Agent('MEK1')
    erk = Agent('ERK2')
    p1 = Phosphorylation(raf, mek,
            evidence=Evidence(text='foo'))
    p2 = Phosphorylation(raf, mek,
            evidence=Evidence(text='bar'))
    p3 = Phosphorylation(raf, mek,
            evidence=Evidence(text='baz'))
    p4 = Phosphorylation(raf, mek,
            evidence=Evidence(text='beep'))
    p5 = Phosphorylation(mek, erk,
            evidence=Evidence(text='foo'))
    p6 = Dephosphorylation(mek, erk,
            evidence=Evidence(text='bar'))
    p7 = Dephosphorylation(mek, erk,
            evidence=Evidence(text='baz'))
    p8 = Dephosphorylation(mek, erk,
            evidence=Evidence(text='beep'))
    p9 = Dephosphorylation(Agent('SRC'), Agent('KRAS'),
                           evidence=Evidence(text='beep'))
    stmts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    pa = Preassembler(hierarchies, stmts=stmts)
    pa.combine_duplicates()
    # The statements come out sorted by their matches_key
    assert(len(pa.unique_stmts) == 4)
    assert(pa.unique_stmts[0].matches(p6)) # MEK dephos ERK
    assert(len(pa.unique_stmts[0].evidence) == 3)
    assert(pa.unique_stmts[1].matches(p9)) # SRC dephos KRAS
    assert(len(pa.unique_stmts[1].evidence) == 1)
    assert(pa.unique_stmts[2].matches(p5)) # MEK phos ERK
    assert(len(pa.unique_stmts[2].evidence) == 1)
    assert(pa.unique_stmts[3].matches(p1)) # RAF phos MEK
    assert(len(pa.unique_stmts[3].evidence) == 4)

def test_superfamily_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'BE': 'RAS'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    st1 = Phosphorylation(src, ras, 'tyrosine', '32')
    st2 = Phosphorylation(src, nras, 'tyrosine', '32')
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # The top-level list should contain only one statement, the gene-level
    # one, supported by the family one.
    assert(len(stmts) == 1)
    assert (stmts[0].equals(st2))
    assert (len(stmts[0].supported_by) == 1)
    assert (stmts[0].supported_by[0].equals(st1))

def test_modification_refinement():
    """A more specific modification statement should be supported by a more
    generic modification statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    st1 = Phosphorylation(src, nras, 'tyrosine', '32')
    st2 = Phosphorylation(src, nras)
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # The top-level list should contain only one statement, the more specific
    # modification, supported by the less-specific modification.
    assert(len(stmts) == 1)
    assert (stmts[0].equals(st1))
    assert (len(stmts[0].supported_by) == 1)
    assert (stmts[0].supported_by[0].equals(st2))

def test_modification_refinement_residue_noenz():
    erbb3 = Agent('Erbb3')
    st1 = Phosphorylation(None, erbb3)
    st2 = Phosphorylation(None, erbb3, 'Y')
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    pa.combine_related()
    assert(len(pa.related_stmts) == 1)

def test_modification_refinement_noenz():
    """A more specific modification statement should be supported by a more
    generic modification statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    st1 = Phosphorylation(src, nras, 'tyrosine', '32')
    st2 = Phosphorylation(None, nras, 'tyrosine', '32')
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # The top-level list should contain only one statement, the more specific
    # modification, supported by the less-specific modification.
    assert(len(stmts) == 1)
    assert (stmts[0].equals(st1))
    assert (stmts[0].equals(st1))
    assert (len(stmts[0].supported_by) == 1)
    assert (stmts[0].supported_by[0].equals(st2))

def test_modification_norefinement_noenz():
    """A more specific modification statement should be supported by a more
    generic modification statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    st1 = Phosphorylation(src, nras)
    st2 = Phosphorylation(None, nras, 'Y', '32',
                          evidence=[Evidence(text='foo')])
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # Modification is less specific, enzyme more specific in st1, therefore
    # these statements shouldn't be combined. 
    assert(len(stmts) == 2)
    assert(len(stmts[1].evidence)==1)

def test_modification_norefinement_subsfamily():
    """A more specific modification statement should be supported by a more
    generic modification statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    ras = Agent('RAS', db_refs = {'BE': 'RAS'})
    st1 = Phosphorylation(src, nras)
    st2 = Phosphorylation(src, ras, 'Y', '32',
                          evidence=[Evidence(text='foo')])
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # Modification is less specific, enzyme more specific in st1, therefore
    # these statements shouldn't be combined. 
    assert(len(stmts) == 2)
    assert(len(stmts[1].evidence)==1)

def test_modification_norefinement_enzfamily():
    """A more specific modification statement should be supported by a more
    generic modification statement."""
    mek = Agent('MEK')
    raf = Agent('RAF')
    braf = Agent('BRAF')
    st1 = Phosphorylation(raf, mek, 'Y', '32',
                          evidence=[Evidence(text='foo')])
    st2 = Phosphorylation(braf, mek)
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # Modification is less specific, enzyme more specific in st1, therefore
    # these statements shouldn't be combined. 
    assert(len(stmts) == 2)
    assert(len(stmts[1].evidence)==1)

def test_bound_condition_refinement():
    """A statement with more specific bound context should be supported by a
    less specific statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    gtp = Agent('GTP', db_refs = {'CHEBI': '15996'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nrasgtp = Agent('NRAS', db_refs = {'HGNC': '7989'},
        bound_conditions=[BoundCondition(gtp, True)])
    st1 = Phosphorylation(src, nras, 'tyrosine', '32')
    st2 = Phosphorylation(src, nrasgtp, 'tyrosine', '32')
    # The top-level list should contain only one statement, the more specific
    # modification, supported by the less-specific modification.
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    assert(len(stmts) == 1)
    assert (stmts[0].equals(st2))
    assert (len(stmts[0].supported_by) == 1)
    assert (stmts[0].supported_by[0].equals(st1))

def test_bound_condition_norefinement():
    """A statement with more specific bound context should be supported by a
    less specific statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    gtp = Agent('GTP', db_refs = {'CHEBI': '15996'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nrasgtp = Agent('NRAS', db_refs = {'HGNC': '7989'},
        bound_conditions=[BoundCondition(gtp, True)])
    st1 = Phosphorylation(src, nras, 'tyrosine', '32')
    st2 = Phosphorylation(src, nrasgtp)
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # The bound condition is more specific in st2 but the modification is less
    # specific. Therefore these statements should not be combined.
    assert(len(stmts) == 2)

def test_complex_refinement():
    ras = Agent('RAS')
    raf = Agent('RAF')
    mek = Agent('MEK')
    st1 = Complex([ras, raf])
    st2 = Complex([mek, ras, raf])
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    pa.combine_related()
    assert(len(pa.unique_stmts) == 2)
    assert(len(pa.related_stmts) == 2)

def test_mod_sites_refinement():
    """A statement with more specific modification context should be supported
    by a less-specific statement."""
    # TODO
    assert True

def test_binding_site_refinement():
    """A statement with information about a binding site for an interaction
    between two proteins should be supported by a statement without this
    information."""
    # TODO
    assert True

def test_activating_substitution_refinement():
    """Should only be refinement if entities are a refinement and all
    fields match."""
    mc1 = MutCondition('12', 'G', 'D')
    mc2 = MutCondition('61', 'Q', 'L')
    nras1 = Agent('NRAS', mutations=[mc1], db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', mutations=[mc2], db_refs = {'HGNC': '7989'})
    ras = Agent('RAS', mutations=[mc1], db_refs={'BE': 'RAS'})
    st1 = ActiveForm(ras, 'gtpbound1', True,
                     evidence=Evidence(text='bar'))
    st2 = ActiveForm(nras1, 'gtpbound1', True,
                     evidence=Evidence(text='foo'))
    st3 = ActiveForm(nras2, 'gtpbound1', True,
                     evidence=Evidence(text='bar'))
    st4 = ActiveForm(nras1, 'gtpbound2', True,
                     evidence=Evidence(text='bar'))
    st5 = ActiveForm(nras1, 'gtpbound1', False,
                     evidence=Evidence(text='bar'))
    assert(st2.refinement_of(st1, hierarchies))
    assert(not st3.refinement_of(st1, hierarchies))
    assert(not st4.refinement_of(st1, hierarchies))
    assert(not st5.refinement_of(st1, hierarchies))

    assert(not st1.refinement_of(st2, hierarchies))
    assert(not st3.refinement_of(st2, hierarchies))
    assert(not st4.refinement_of(st2, hierarchies))
    assert(not st5.refinement_of(st2, hierarchies))

    assert(not st1.refinement_of(st3, hierarchies))
    assert(not st2.refinement_of(st3, hierarchies))
    assert(not st4.refinement_of(st3, hierarchies))
    assert(not st5.refinement_of(st3, hierarchies))

    assert(not st1.refinement_of(st4, hierarchies))
    assert(not st2.refinement_of(st4, hierarchies))
    assert(not st3.refinement_of(st4, hierarchies))
    assert(not st5.refinement_of(st4, hierarchies))

    assert(not st1.refinement_of(st5, hierarchies))
    assert(not st2.refinement_of(st5, hierarchies))
    assert(not st3.refinement_of(st5, hierarchies))
    assert(not st4.refinement_of(st5, hierarchies))

def test_translocation():
    st1 = Translocation(Agent('AKT'), None, None)
    st2 = Translocation(Agent('AKT'), None, 'plasma membrane')
    st3 = Translocation(Agent('AKT'), None, 'nucleus')
    pa = Preassembler(hierarchies, stmts=[st1, st2, st3])
    pa.combine_related()
    assert(len(pa.related_stmts) == 2)

def test_render_stmt_graph():
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    mek = Agent('MEK', db_refs={'BE':'MEK'})
    # Statements
    p0 = Phosphorylation(braf, mek)
    p1 = Phosphorylation(braf, mek1)
    p2 = Phosphorylation(braf, mek1, position='218')
    p3 = Phosphorylation(braf, mek1, position='222')
    p4 = Phosphorylation(braf, mek1, 'serine')
    p5 = Phosphorylation(braf, mek1, 'serine', '218')
    p6 = Phosphorylation(braf, mek1, 'serine', '222')
    stmts = [p0, p1, p2, p3, p4, p5, p6]
    pa = Preassembler(hierarchies, stmts=stmts)
    pa.combine_related()
    graph = render_stmt_graph(pa.related_stmts)
    # One node for each statement
    assert len(graph.nodes()) == 7
    # Edges:
    # p0 supports p1-p6 = 6 edges
    # p1 supports p2-p6 = 5 edges
    # p2 supports p5 = 1 edge
    # p3 supports p6 = 1 edge
    # p4 supports p5-p6 = 2 edges
    # (p5 and p6 support none--they are top-level)
    # 6 + 5 + 1 + 1 + 2 = 15 edges
    assert len(graph.edges()) == 15

def test_flatten_evidence_hierarchy():
    braf = Agent('BRAF')
    mek = Agent('MAP2K1')
    st1 = Phosphorylation(braf, mek, evidence=[Evidence(text='foo')])
    st2 = Phosphorylation(braf, mek, 'S', '218',
                          evidence=[Evidence(text='bar')])
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    pa.combine_related()
    assert len(pa.related_stmts) == 1
    flattened = flatten_evidence(pa.related_stmts)
    assert len(flattened) == 1
    top_stmt = flattened[0]
    assert len(top_stmt.evidence) == 2
    assert 'bar' in [e.text for e in top_stmt.evidence]
    assert 'foo' in [e.text for e in top_stmt.evidence]
    assert len(top_stmt.supported_by) == 1
    supporting_stmt = top_stmt.supported_by[0]
    assert len(supporting_stmt.evidence) == 1
    assert supporting_stmt.evidence[0].text == 'foo'

def test_flatten_stmts():
    st1 = Phosphorylation(Agent('MAP3K5'), Agent('RAF1'), 'S', '338')
    st2 = Phosphorylation(None, Agent('RAF1'), 'S', '338')
    st3 = Phosphorylation(None, Agent('RAF1'))
    st4 = Phosphorylation(Agent('PAK1'), Agent('RAF1'), 'S', '338')
    st5 = Phosphorylation(None, Agent('RAF1'), evidence=Evidence(text='foo'))
    pa = Preassembler(hierarchies, stmts=[st1, st2, st3, st4, st5])
    pa.combine_duplicates()
    pa.combine_related()
    assert(len(pa.related_stmts) == 2)
    assert(len(flatten_stmts(pa.unique_stmts)) == 4)
    assert(len(flatten_stmts(pa.related_stmts)) == 4)

def test_complex_refinement_order():
    st1 = Complex([Agent('MED23'), Agent('ELK1')])
    st2 = Complex([Agent('ELK1', mods=[ModCondition('phosphorylation')]),
                   Agent('MED23')])
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    pa.combine_duplicates()
    pa.combine_related()
    assert(len(pa.related_stmts) == 1)

if __name__ == '__main__':
    test_modification_refinement_residue_noenz()


