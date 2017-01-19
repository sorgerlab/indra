from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
import itertools
from indra.util import unicode_strs
from indra.tools import expand_families as ef
from indra.preassembler.hierarchy_manager import hierarchies
from indra.statements import Agent, Phosphorylation, Complex, Activation

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
    expanded_stmts = exp.expand_families([st])
    assert len(expanded_stmts) == 2
    # Test for case involving None for one of the agents
    st = Phosphorylation(None, akt)
    expanded_stmts = exp.expand_families([st])
    assert len(expanded_stmts) == 3
    # Statement with two families: 3 Rafs x 2 Meks
    st = Phosphorylation(raf, mek, 'S', '202')
    expanded_stmts = exp.expand_families([st])
    assert len(expanded_stmts) == 6
    # Test also for case involving both family and complex relationships
    st = Phosphorylation(ampk, mek)
    expanded_stmts = exp.expand_families([st])
    assert len(expanded_stmts) == 14

def test_complexes_from_hierarchy():
    exp = ef.Expander(hierarchies)
    complexes = exp.complexes_from_hierarchy()
    keys = [c.matches_key() for c in complexes]
    probe_stmt = Complex([Agent('AMPK_alpha', db_refs={'BE':'AMPK_alpha'}),
                          Agent('AMPK_beta', db_refs={'BE':'AMPK_beta'}),
                          Agent('AMPK_gamma', db_refs={'BE':'AMPK_gamma'})])
    assert probe_stmt.matches_key() in keys

def test_expanded_complexes_from_hierarchy():
    exp = ef.Expander(hierarchies)
    complexes = exp.expanded_complexes_from_hierarchy()
    stmt_ag_names = []
    for stmt in complexes:
        sorted_names = tuple(sorted([ag.name for ag in stmt.agent_list()]))
        stmt_ag_names.append(sorted_names)
    ampk_alphas = ('PRKAA1', 'PRKAA2')
    ampk_betas = ('PRKAB1', 'PRKAB2')
    ampk_gammas = ('PRKAG1', 'PRKAG2', 'PRKAG3')
    for alpha, beta, gamma in itertools.product(ampk_alphas, ampk_betas,
                                                ampk_gammas):
        assert tuple(sorted((alpha, beta, gamma))) in stmt_ag_names

def test_db_ref_keys():
    # test that expanded families get TEXT, UP, HGNC keys in their db_refs
    # Get the Expander
    exp = ef.Expander(hierarchies)
    # Declare some agents
    grb2 = Agent('GRB2',
                 db_refs={'TEXT': 'Grb2', 'UP': 'P62993', 'HGNC': '4566'})
    shc = Agent('SHC', db_refs={'BE':'SHC'})
    # Test case where one agent is a family and the other is a gene
    st = Activation(grb2, shc)
    expanded_stmts = exp.expand_families([st])
    for st in expanded_stmts:
        for agent in st.agent_list():
            if agent is not None:
                assert set(list(agent.db_refs)) == \
                set(['TEXT','UP','HGNC'])
    # Test for case involving None for one of the agents
    st = Phosphorylation(None, shc)
    expanded_stmts = exp.expand_families([st])
    for st in expanded_stmts:
        for agent in st.agent_list():
            if agent is not None:
                assert set(list(agent.db_refs)) == \
                set(['TEXT','UP','HGNC'])
    # Statement with two families: 4x4 SHC
    st = Activation(shc, shc)
    expanded_stmts = exp.expand_families([st])
    for st in expanded_stmts:
        for agent in st.agent_list():
            if agent is not None:
                assert set(list(agent.db_refs)) == \
                set(['TEXT','UP','HGNC'])

def test_get_children():
    exp = ef.Expander(hierarchies)
    raf = Agent('RAF', db_refs={'BE':'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC':'1097'})
    # Look up RAF
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
    assert len(mapks) == 12
    assert ('HGNC', 'MAPK1') in mapks
    assert ('HGNC', 'MAPK9') in mapks
    assert ('BE', 'ERK') in mapks
    assert ('BE', 'JNK') in mapks
    assert unicode_strs(mapks)
    # Now do the same expansion with a namespace filter
    mapks = exp.get_children(mapk, ns_filter='HGNC')
    assert unicode_strs(mapks)
    assert len(mapks) == 9
    assert ('HGNC', 'MAPK3') in mapks
    assert ('HGNC', 'MAPK10') in mapks
    assert ('BE', 'ERK') not in mapks
    # Make sure we can also do this in a case involving both family and complex
    # relationships
    ampk = Agent('AMPK', db_refs={'BE':'AMPK'})
    ampks = exp.get_children(ampk, ns_filter=None)
    assert len(ampks) == 22
    ampks = exp.get_children(ampk, ns_filter='HGNC')
    assert len(ampks) == 7
    # Test that the default filter is HGNC
    ampks = exp.get_children(ampk)
    assert len(ampks) == 7
    ag_none = None
    none_children = exp.get_children(ag_none)
    assert isinstance(none_children, list)
    assert len(none_children) == 0
