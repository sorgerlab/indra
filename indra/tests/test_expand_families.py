from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
from indra.preassembler.hierarchy_manager import hierarchies

from indra.statements import Agent, Phosphorylation
from indra.tools import expand_families as ef
from indra.util import unicode_strs

def test_component_children_lookup():
    raf = Agent('RAF', db_refs={'BE':'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC':'1097'})
    # Look up RAF
    exp = ef.Expander(hierarchies)
    rafs = exp.get_children(raf)
    # Should get three family members
    assert isinstance(rafs, list)
    assert len(rafs) == 3
    assert unicode_strs(rafs)
    # The lookup of a gene-level entity should not return any additional
    # entities
    brafs = exp.get_children(braf)
    assert isinstance(brafs, list)
    assert len(brafs) == 0
    assert unicode_strs(brafs)
    # The lookup for a top-level family (e.g., MAPK, which has as children
    # both the intermediate family ERK as well as all MAPK1-15 members)
    # should not return the intermediate families when a filter is applied.
    mapk = Agent('MAPK', db_refs={'BE':'MAPK'})
    mapks = exp.get_children(mapk, ns_filter=None)
    assert len(mapks) == 3
    assert ('HGNC', 'MAPK1') in mapks
    assert ('HGNC', 'MAPK3') in mapks
    assert ('BE', 'ERK') in mapks
    assert unicode_strs(mapks)
    # Now do the same expansion with a namespace filter
    mapks = exp.get_children(mapk, ns_filter='HGNC')
    assert unicode_strs(mapks)
    assert len(mapks) == 2
    assert ('HGNC', 'MAPK1') in mapks
    assert ('HGNC', 'MAPK3') in mapks
    # Make sure we can also do this in a case involving both family and complex
    # relationships
    ampk = Agent('AMPK', db_refs={'BE':'AMPK'})
    ampks = exp.get_children(ampk, ns_filter=None)
    assert len(ampks) == 10
    ampks = exp.get_children(ampk, ns_filter='HGNC')
    assert len(ampks) == 7
    # Test that the default filter is HGNC
    ampks = exp.get_children(ampk)
    assert len(ampks) == 7
    ag_none = None
    none_children = exp.get_children(ag_none)
    assert isinstance(none_children, list)
    assert len(none_children) == 0

def test_expand_families():
    # Get the Expander
    exp = ef.Expander(hierarchies)
    # Declare some agents
    akt = Agent('AKT', db_refs={'BE':'AKT'})
    raf = Agent('RAF', db_refs={'BE':'RAF'})
    mek = Agent('MEK', db_refs={'BE':'MEK'})
    mapk1 = Agent('MAPK1', db_refs={'BE':'MAPK1'})
    ampk = Agent('AMPK', db_refs={'BE':'AMPK'})
    # Test case where one agent is a family and the other is a gene
    st = Phosphorylation(mek, mapk1)
    import ipdb; ipdb.set_trace()
    expanded_stmts = exp.expand_families([st])
    assert len(expanded_stmts) == 2
    # Test for case involving None for one of the agents
    st = Phosphorylation(None, akt)
    expanded_stmts = exp.expand_families([st])
    assert len(expanded_stmts) == 3

    st = Phosphorylation(raf, mek, 'S', '202')
    expanded_stmts = exp.expand_families([st])
    # 3 Rafs x 2 Meks
    assert len(expanded_stmts) == 6
    # Test also for case involving both family and complex relationships
    st = Phosphorylation(ampk, mek)
    expanded_stmts = exp.expand_families([st])
    assert len(expanded_stmts) == 14

if __name__ == '__main__':
    test_expand_families()
    #test_component_children_lookup()
