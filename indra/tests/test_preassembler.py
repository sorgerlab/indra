import os
from indra.preassembler import Preassembler
from indra.trips import trips_api
from indra.statements import Agent, Phosphorylation, BoundCondition, \
                             Dephosphorylation, Evidence
from indra.preassembler.hierarchy_manager import HierarchyManager

entity_file_path = os.path.join(os.path.dirname(__file__),
                    '..', 'preassembler', 'entity_hierarchy.rdf')
mod_file_path = os.path.join(os.path.dirname(__file__),
                    '..', 'preassembler', 'modification_hierarchy.rdf')
eh = HierarchyManager(entity_file_path)
mh = HierarchyManager(mod_file_path)

"""
def test_from_text():
    sentences = ['Src phosphorylates Ras, bound to GTP, at Tyr32.', 
                'Src phosphorylates NRAS at Tyr32.',
                'Src phosphorylates NRAS that is bound to GTP.']
    pa = Preassembler()
    for s in sentences:
        tp = trips_api.process_text(s)
        pa.add_statements(tp.statements)
"""

def test_duplicates():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FA': '03663'})
    st1 = Phosphorylation(src, ras, 'Phosphorylation', None)
    st2 = Phosphorylation(src, ras, 'Phosphorylation', None)
    pa = Preassembler(eh, mh, [st1, st2])
    pa.combine_duplicates()
    assert(len(pa.unique_stmts) == 1)

def test_combine_duplicates():
    raf = Agent('RAF1')
    mek = Agent('MEK1')
    erk = Agent('ERK2')
    p1 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='foo'))
    p2 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='bar'))
    p3 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='baz'))
    p4 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='beep'))
    p5 = Phosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='foo'))
    p6 = Dephosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='bar'))
    p7 = Dephosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='baz'))
    p8 = Dephosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='beep'))
    p9 = Dephosphorylation(Agent('SRC'), Agent('KRAS'),
                         'Phosphorylation', None, evidence=Evidence(text='beep'))
    stmts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    pa = Preassembler(eh, mh, stmts)
    pa.combine_duplicates()
    # The statements come out sorted by their matches_key
    assert(len(pa.unique_stmts) == 4)
    assert(pa.unique_stmts[0] == p5) # MEK phos ERK
    assert(len(pa.unique_stmts[0].evidence) == 1)
    assert(pa.unique_stmts[1] == p1) # RAF phos MEK
    assert(len(pa.unique_stmts[1].evidence) == 4)
    assert(pa.unique_stmts[2] == p6) # MEK dephos ERK
    assert(len(pa.unique_stmts[2].evidence) == 3)
    assert(pa.unique_stmts[3] == p9) # SRC dephos KRAS
    assert(len(pa.unique_stmts[3].evidence) == 1)

def test_superfamily_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FA': '03663'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    st1 = Phosphorylation(src, ras, 'PhosphorylationTyrosine', '32')
    st2 = Phosphorylation(src, nras, 'PhosphorylationTyrosine', '32')
    pa = Preassembler(eh, mh, [st1, st2])
    stmts = pa.combine_related()
    # The top-level list should contain only one statement, the gene-level
    # one, supported by the family one.
    assert(len(stmts) == 1)
    assert (stmts[0] == st2)
    assert (stmts[0].supported_by == [st1])


def test_modification_refinement():
    """A more specific modification statement should be supported by a more
    generic modification statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    st1 = Phosphorylation(src, nras, 'PhosphorylationTyrosine', '32')
    st2 = Phosphorylation(src, nras, 'Phosphorylation', None)
    pa = Preassembler(eh, mh, [st1, st2])
    stmts = pa.combine_related()
    # The top-level list should contain only one statement, the more specific
    # modification, supported by the less-specific modification.
    assert(len(stmts) == 1)
    assert (stmts[0] == st1)
    assert (stmts[0].supported_by == [st2])

def test_bound_condition_refinement():
    """A statement with more specific bound context should be supported by a
    less specific statement."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    gtp = Agent('GTP', db_refs = {'CHEBI': '15996'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nrasgtp = Agent('NRAS', db_refs = {'HGNC': '7989'},
        bound_conditions=[BoundCondition(gtp, True)])
    st1 = Phosphorylation(src, nras, 'PhosphorylationTyrosine', '32')
    st2 = Phosphorylation(src, nrasgtp, 'PhosphorylationTyrosine', '32')
    # The top-level list should contain only one statement, the more specific
    # modification, supported by the less-specific modification.
    pa = Preassembler(eh, mh, [st1, st2])
    stmts = pa.combine_related()
    assert(len(stmts) == 1)
    assert (stmts[0] == st2)
    assert (stmts[0].supported_by == [st1])

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

