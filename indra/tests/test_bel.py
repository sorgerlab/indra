from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra import bel
from indra.util import unicode_strs
from rdflib.term import URIRef
from indra.bel.processor import BelProcessor

concept_prefix = 'http://www.openbel.org/bel/namespace//'
entity_prefix = 'http://www.openbel.org/bel/'

def test_bel_ndex_query():
    bp = bel.process_ndex_neighborhood(['NFKB1'])
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

