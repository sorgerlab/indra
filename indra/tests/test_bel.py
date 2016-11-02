from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra import bel
from indra.util import unicode_strs
from rdflib.term import URIRef
from indra.bel.processor import BelProcessor

concept_prefix = 'http://www.openbel.org/bel/namespace//'
entity_prefix = 'http://www.openbel.org/bel/'

path_this = os.path.dirname(os.path.abspath(__file__))
test_rdf = os.path.join(path_this, 'bel_rdfs', 'NFKB1_neighborhood.rdf')
with open(test_rdf, 'rb') as fh:
    rdf_str = fh.read()

def assert_pmids(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid is not None:
                assert(ev.pmid.isdigit())

def test_bel_ndex_query():
    bp = bel.process_ndex_neighborhood(['NFKB1'])
    assert_pmids(bp.statements)
    unicode_strs(bp.statements)

def test_process_belrdf():
    bp = bel.process_belrdf(rdf_str)
    assert_pmids(bp.statements)
    unicode_strs(bp.statements)

def test_get_agent_up_from_hgnc():
    hgnc_sym = 'MAPK1'
    concept = concept_prefix + hgnc_sym
    entity = entity_prefix + 'p_HGNC_' + hgnc_sym
    ag = BelProcessor.get_agent(concept, entity)
    assert ag.name == 'MAPK1'
    assert ag.db_refs.get('HGNC') == '6871'
    assert ag.db_refs.get('UP') == 'P28482'
    assert unicode_strs((concept, entity, ag))

def test_get_agent_hgnc_up_from_egid():
    entrez_id = '5594'
    concept = concept_prefix + entrez_id
    entity = entity_prefix + 'p_EGID_' + entrez_id
    ag = BelProcessor.get_agent(concept, entity)
    assert ag.name == 'MAPK1'
    assert ag.db_refs.get('EGID') == entrez_id
    assert ag.db_refs.get('HGNC') == '6871'
    assert ag.db_refs.get('UP') == 'P28482'
    assert unicode_strs((concept, entity, ag))

