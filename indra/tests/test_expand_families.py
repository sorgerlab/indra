import itertools
from indra.util import unicode_strs
from indra.tools import expand_families as ef
from indra.ontology.bio import bio_ontology
from indra.statements import Agent, Phosphorylation, Complex, Activation


def test_expand_families():
    # Get the Expander
    exp = ef.Expander(bio_ontology)
    # Declare some agents
    akt = Agent('AKT', db_refs={'FPLX': 'AKT'})
    raf = Agent('RAF', db_refs={'FPLX': 'RAF'})
    mek = Agent('MEK', db_refs={'FPLX': 'MEK'})
    mapk1 = Agent('MAPK1', db_refs={'FPLX': 'MAPK1'})
    ampk = Agent('AMPK', db_refs={'FPLX': 'AMPK'})
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
    exp = ef.Expander(bio_ontology)
    complexes = exp.complexes_from_hierarchy()
    keys = [c.matches_key() for c in complexes]
    probe_stmt = Complex([Agent('AMPK_alpha', db_refs={'FPLX': 'AMPK_alpha'}),
                          Agent('AMPK_beta', db_refs={'FPLX': 'AMPK_beta'}),
                          Agent('AMPK_gamma', db_refs={'FPLX': 'AMPK_gamma'})])
    assert probe_stmt.matches_key() in keys


def test_expanded_complexes_from_hierarchy():
    exp = ef.Expander(bio_ontology)
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
    exp = ef.Expander(bio_ontology)
    # Declare some agents
    grb2 = Agent('GRB2',
                 db_refs={'TEXT': 'Grb2', 'UP': 'P62993', 'HGNC': '4566'})
    shc = Agent('SHC', db_refs={'FPLX': 'SHC'})
    # Test case where one agent is a family and the other is a gene
    st = Activation(grb2, shc)
    expanded_stmts = exp.expand_families([st])
    for st in expanded_stmts:
        for agent in st.agent_list():
            if agent is not None:
                assert set(list(agent.db_refs)) == {'TEXT', 'UP', 'HGNC'}
    # Test for case involving None for one of the agents
    st = Phosphorylation(None, shc)
    expanded_stmts = exp.expand_families([st])
    for st in expanded_stmts:
        for agent in st.agent_list():
            if agent is not None:
                assert set(list(agent.db_refs)) == {'TEXT', 'UP', 'HGNC'}
    # Statement with two families: 4x4 SHC
    st = Activation(shc, shc)
    expanded_stmts = exp.expand_families([st])
    for st in expanded_stmts:
        for agent in st.agent_list():
            if agent is not None:
                assert set(list(agent.db_refs)) == {'TEXT', 'UP', 'HGNC'}
