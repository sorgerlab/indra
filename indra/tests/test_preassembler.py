from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra.preassembler import Preassembler, render_stmt_graph, \
                               flatten_evidence, flatten_stmts
from indra.preassembler.hierarchy_manager import HierarchyManager
from indra.sources import trips
from indra.statements import Agent, Phosphorylation, BoundCondition, \
                             Dephosphorylation, Evidence, ModCondition, \
                             ActiveForm, MutCondition, Complex, \
                             Translocation, Activation, Inhibition, \
                             Deacetylation, Conversion, Concept, Influence, \
                             IncreaseAmount, DecreaseAmount
from indra.preassembler.hierarchy_manager import hierarchies

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
            evidence=Evidence(text='foo2'))
    p6 = Dephosphorylation(mek, erk,
            evidence=Evidence(text='bar2'))
    p7 = Dephosphorylation(mek, erk,
            evidence=Evidence(text='baz2'))
    p8 = Dephosphorylation(mek, erk,
            evidence=Evidence(text='beep2'))
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

def test_combine_evidence_exact_duplicates():
    raf = Agent('RAF1')
    mek = Agent('MEK1')
    p1 = Phosphorylation(raf, mek,
            evidence=Evidence(text='foo'))
    p2 = Phosphorylation(raf, mek,
            evidence=Evidence(text='bar'))
    p3 = Phosphorylation(raf, mek,
            evidence=Evidence(text='bar'))
    stmts = [p1, p2, p3]
    pa = Preassembler(hierarchies, stmts=stmts)
    pa.combine_duplicates()
    # The statements come out sorted by their matches_key
    assert(len(pa.unique_stmts) == 1)
    assert(len(pa.unique_stmts[0].evidence) == 2)
    assert(set(ev.text for ev in pa.unique_stmts[0].evidence) ==
           set(['foo', 'bar']))

def test_superfamily_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FPLX': 'RAS'})
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


def test_superfamily_refinement_isa_or_partof():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    prkag1 = Agent('PRKAG1', db_refs = {'HGNC': '9385'})
    ampk = Agent('AMPK', db_refs = {'FPLX': 'AMPK'})
    st1 = Phosphorylation(src, ampk, 'tyrosine', '32')
    st2 = Phosphorylation(src, prkag1, 'tyrosine', '32')
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
    assert (len(stmts[0].supported_by) == 1)
    assert (stmts[0].supported_by[0].equals(st2))
    assert (stmts[0].supported_by[0].supports[0].equals(st1))

def test_modification_refinement_noenz2():
    """A more specific modification statement should be supported by a more
    generic modification statement.

    Similar to test_modification_refinement_noenz for statements where one
    argument is associated with a component in the hierarchy (SIRT1 in this
    case) but the other is not (BECN1).
    """
    sirt1 = Agent('SIRT1', db_refs={'HGNC':'14929', 'UP':'Q96EB6',
                                    'TEXT':'SIRT1'})
    becn1 = Agent('BECN1', db_refs={'HGNC': '1034', 'UP': 'Q14457',
                                    'TEXT': 'Beclin 1'})
    st1 = Deacetylation(sirt1, becn1)
    st2 = Deacetylation(None, becn1)
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    # The top-level list should contain only one statement, the more specific
    # modification, supported by the less-specific modification.
    assert (len(stmts) == 1)
    assert (stmts[0].equals(st1))
    assert (len(stmts[0].supported_by) == 1)
    assert (stmts[0].supported_by[0].equals(st2))
    assert (stmts[0].supported_by[0].supports[0].equals(st1))

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
    ras = Agent('RAS', db_refs = {'FPLX': 'RAS'})
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

def test_bound_condition_deep_refinement():
    """A statement with more specific bound context should be supported by a
    less specific statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    gtp1 = Agent('GTP', db_refs = {'CHEBI': '15996'})
    gtp2 = Agent('GTP', mods=[ModCondition('phosphorylation')],
                 db_refs = {'CHEBI': '15996'})
    nrasgtp1 = Agent('NRAS', db_refs = {'HGNC': '7989'},
        bound_conditions=[BoundCondition(gtp1, True)])
    nrasgtp2 = Agent('NRAS', db_refs = {'HGNC': '7989'},
        bound_conditions=[BoundCondition(gtp2, True)])
    st1 = Phosphorylation(src, nrasgtp1, 'tyrosine', '32')
    st2 = Phosphorylation(src, nrasgtp2, 'tyrosine', '32')
    # The top-level list should contain only one statement, the more specific
    # modification, supported by the less-specific modification.
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related()
    assert(len(stmts) == 1)
    assert (stmts[0].equals(st2))
    assert (len(stmts[0].supported_by) == 1)
    assert (stmts[0].supported_by[0].equals(st1))

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

def test_complex_agent_refinement():
    ras = Agent('RAS')
    raf1 = Agent('RAF', mods=[ModCondition('ubiquitination', None, None, True)])
    raf2 = Agent('RAF', mods=[ModCondition('ubiquitination', None, None, False)])
    st1 = Complex([ras, raf1])
    st2 = Complex([ras, raf2])
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
    ras = Agent('RAS', mutations=[mc1], db_refs={'FPLX': 'RAS'})
    st1 = ActiveForm(ras, 'gtpbound', True,
                     evidence=Evidence(text='bar'))
    st2 = ActiveForm(nras1, 'gtpbound', True,
                     evidence=Evidence(text='foo'))
    st3 = ActiveForm(nras2, 'gtpbound', True,
                     evidence=Evidence(text='bar'))
    st4 = ActiveForm(nras1, 'phosphatase', True,
                     evidence=Evidence(text='bar'))
    st5 = ActiveForm(nras1, 'gtpbound', False,
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

def test_grounding_aggregation():
    braf1 = Agent('BRAF', db_refs={'TEXT': 'braf', 'HGNC': '1097'})
    braf2 = Agent('BRAF', db_refs={'TEXT': 'BRAF'})
    braf3 = Agent('BRAF', db_refs={'TEXT': 'Braf', 'UP': 'P15056'})
    st1 = Phosphorylation(None, braf1)
    st2 = Phosphorylation(None, braf2)
    st3 = Phosphorylation(None, braf3)
    pa = Preassembler(hierarchies, stmts=[st1, st2, st3])
    unique_stmts = pa.combine_duplicates()
    assert(len(unique_stmts) == 3)

def test_grounding_aggregation_complex():
    mek = Agent('MEK')
    braf1 = Agent('BRAF', db_refs={'TEXT': 'braf', 'HGNC': '1097'})
    braf2 = Agent('BRAF', db_refs={'TEXT': 'BRAF', 'dummy': 'dummy'})
    braf3 = Agent('BRAF', db_refs={'TEXT': 'Braf', 'UP': 'P15056'})
    st1 = Complex([mek, braf1])
    st2 = Complex([braf2, mek])
    st3 = Complex([mek, braf3])
    pa = Preassembler(hierarchies, stmts=[st1, st2, st3])
    unique_stmts = pa.combine_duplicates()
    assert(len(unique_stmts) == 3)

def test_render_stmt_graph():
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    mek = Agent('MEK', db_refs={'FPLX':'MEK'})
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
    graph = render_stmt_graph(pa.related_stmts, reduce=False)
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

def test_activation_refinement():
    subj = Agent('alcohol', db_refs={'CHEBI': 'CHEBI:16236',
                                     'HMDB': 'HMDB00108',
                                     'PUBCHEM': '702',
                                     'TEXT': 'alcohol'})
    obj = Agent('endotoxin', db_refs={'TEXT': 'endotoxin'})
    st1 = Inhibition(subj, obj)
    st2 = Activation(subj, obj)
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    pa.combine_duplicates()
    assert(len(pa.unique_stmts) == 2)
    pa.combine_related()
    assert(len(pa.related_stmts) == 2)

def test_homodimer_refinement():
    egfr = Agent('EGFR')
    erbb = Agent('ERBB2')
    st1 = Complex([erbb, erbb])
    st2 = Complex([erbb, egfr])
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    pa.combine_duplicates()
    assert(len(pa.unique_stmts) == 2)
    pa.combine_related()
    assert(len(pa.related_stmts) == 2)

def test_return_toplevel():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    st1 = Phosphorylation(src, nras, 'tyrosine', '32')
    st2 = Phosphorylation(src, nras)
    pa = Preassembler(hierarchies, stmts=[st1, st2])
    stmts = pa.combine_related(return_toplevel=True)
    assert(len(stmts) == 1)
    assert(len(stmts[0].supported_by) == 1)
    assert(len(stmts[0].supported_by[0].supports) == 1)
    stmts = pa.combine_related(return_toplevel=False)
    assert(len(stmts) == 2)
    ix = 1 if stmts[0].residue else 0
    assert(len(stmts[1-ix].supported_by) == 1)
    assert(len(stmts[1-ix].supported_by[0].supports) == 1)
    assert(len(stmts[ix].supports) == 1)
    assert(len(stmts[ix].supports[0].supported_by) == 1)

def test_multiprocessing():
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    mek = Agent('MEK', db_refs={'FPLX':'MEK'})
    # Statements
    p0 = Phosphorylation(braf, mek)
    p1 = Phosphorylation(braf, mek1)
    p2 = Phosphorylation(braf, mek1, position='218')
    p3 = Phosphorylation(braf, mek1, position='222')
    p4 = Phosphorylation(braf, mek1, 'serine')
    p5 = Phosphorylation(braf, mek1, 'serine', '218')
    p6 = Phosphorylation(braf, mek1, 'serine', '222')
    p7 = Dephosphorylation(braf, mek1)
    stmts = [p0, p1, p2, p3, p4, p5, p6, p7]
    pa = Preassembler(hierarchies, stmts=stmts)
    # Size cutoff set to a low number so that one group will run remotely
    # and one locally
    toplevel = pa.combine_related(return_toplevel=True, poolsize=1,
                                  size_cutoff=2)
    assert len(toplevel) == 3, 'Got %d toplevel statements.' % len(toplevel)

def test_conversion_refinement():
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    hras = Agent('HRAS', db_refs={'HGNC': '5173'})
    gtp = Agent('GTP')
    gdp = Agent('GDP')
    st1 = Conversion(ras, gtp, gdp)
    st2 = Conversion(hras, gtp, gdp)
    st3 = Conversion(hras, [gtp, gdp], gdp)
    st4 = Conversion(hras, [gdp, gtp], gdp)
    pa = Preassembler(hierarchies, stmts=[st1, st2, st3, st4])
    toplevel_stmts = pa.combine_related()
    assert(len(toplevel_stmts) == 2)


def test_influence_duplicate():
    gov = 'UN/entities/human/government/government_entity'
    agr = 'UN/entities/natural/crop_technology'
    cgov = Concept('government', db_refs={'UN': [(gov, 1.0)]})
    cagr = Concept('agriculture', db_refs={'UN': [(agr, 1.0)]})
    stmt1 = Influence(cgov, cagr, evidence=[Evidence(source_api='eidos1')])
    stmt2 = Influence(cagr, cgov, evidence=[Evidence(source_api='eidos2')])
    stmt3 = Influence(cgov, cagr, evidence=[Evidence(source_api='eidos3')])
    eidos_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/eidos/eidos_ontology.rdf')
    hm = HierarchyManager(eidos_ont, True, True)
    hierarchies = {'entity': hm}
    pa = Preassembler(hierarchies, [stmt1, stmt2, stmt3])
    unique_stmts = pa.combine_duplicates()
    assert len(unique_stmts) == 2
    assert len(unique_stmts[0].evidence) == 2
    assert len(unique_stmts[1].evidence) == 1
    sources = [e.source_api for e in unique_stmts[0].evidence]
    assert set(sources) == set(['eidos1', 'eidos3'])


def test_influence_refinement():
    tran = 'UN/entities/human/infrastructure/transportation'
    truck = 'UN/entities/human/infrastructure/transportation/' + \
        'transportation_methods'
    agr = 'UN/entities/human/livelihood'
    ctran = Concept('transportation', db_refs={'UN': [(tran, 1.0)]})
    ctruck = Concept('trucking', db_refs={'UN': [(truck, 1.0)]})
    cagr = Concept('agriculture', db_refs={'UN': [(agr, 1.0)]})
    stmt1 = Influence(ctran, cagr, evidence=[Evidence(source_api='eidos1')])
    stmt2 = Influence(ctruck, cagr, evidence=[Evidence(source_api='eidos2')])
    stmt3 = Influence(cagr, ctran, evidence=[Evidence(source_api='eidos3')])
    eidos_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/eidos/eidos_ontology.rdf')
    hm = HierarchyManager(eidos_ont, True, True)
    hierarchies = {'entity': hm}
    pa = Preassembler(hierarchies, [stmt1, stmt2, stmt3])
    rel_stmts = pa.combine_related()
    assert len(rel_stmts) == 2
    truck_stmt = [st for st in rel_stmts if st.subj.name == 'trucking'][0]
    assert len(truck_stmt.supported_by) == 1
    assert truck_stmt.supported_by[0].subj.name == 'transportation'


def test_find_contradicts():
    st1 = Inhibition(Agent('a'), Agent('b'))
    st2 = Activation(Agent('a'), Agent('b'))
    st3 = IncreaseAmount(Agent('a'), Agent('b'))
    st4 = DecreaseAmount(Agent('a'), Agent('b'))
    pa = Preassembler(hierarchies, [st1, st2, st3, st4])
    contradicts = pa.find_contradicts()
    assert len(contradicts) == 2
    for s1, s2 in contradicts:
        assert {s1.uuid, s2.uuid} in ({st1.uuid, st2.uuid},
                                      {st3.uuid, st4.uuid})


def test_preassemble_related_complex():
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    kras = Agent('KRAS', db_refs={'HGNC': '6407'})
    hras = Agent('HRAS', db_refs={'HGNC': '5173'})
    st1 = Complex([kras, hras])
    st2 = Complex([kras, ras])
    st3 = Complex([hras, kras])
    st4 = Complex([ras, kras])
    pa = Preassembler(hierarchies, [st1, st2, st3, st4])
    uniq = pa.combine_duplicates()
    assert len(uniq) == 2
    top = pa.combine_related()
    assert len(top) == 1

