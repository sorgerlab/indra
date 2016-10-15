from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
from indra.preassembler.hierarchy_manager import hierarchies

from indra.statements import Agent, Phosphorylation
from indra.tools import expand_families as ef

def test_component_children_lookup():
    raf = Agent('RAF', db_refs={'BE':'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC':'1097'})
    # Look up RAF
    import ipdb; ipdb.set_trace()
    exp = ef.Expander(hierarchies)
    rafs = exp.get_children(raf)
    assert isinstance(rafs, list)
    assert len(rafs) == 3
    # The lookup of a gene-level entity should not return any additional
    # entities
    brafs = exp.get_children(braf)
    assert isinstance(brafs, list)
    assert len(brafs) == 0
    # The lookup for a top-level family (e.g., MAPK, which has as children
    # both the intermediate family ERK as well as all MAPK1-15 members)
    # should not return the intermediate families when a filter is applied.
    #mapk = Agent('MAPK', db_refs={'BE':'MAPK'})
    #mapks = exp.get_children(mapks)

def test_expand_families():
    raf = Agent('RAF', db_refs={'BE':'RAF'})
    mek = Agent('MEK', db_refs={'BE':'MEK'})
    st = Phosphorylation(raf, mek)
    exp = ef.Expander(hierarchies)
    expanded_stmts = exp.expand_families([st])

if __name__ == '__main__':
    test_component_children_lookup()
