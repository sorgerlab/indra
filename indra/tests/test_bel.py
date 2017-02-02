from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra import bel
from indra.util import unicode_strs
from rdflib.term import URIRef
from indra.bel.processor import BelProcessor
from indra.statements import RegulateAmount

concept_prefix = 'http://www.openbel.org/bel/namespace//'
entity_prefix = 'http://www.openbel.org/bel/'

path_this = os.path.dirname(os.path.abspath(__file__))
test_rdf_nfkb = os.path.join(path_this, 'bel_rdfs', 'NFKB1_neighborhood.rdf')
test_rdf_myc = os.path.join(path_this, 'bel_rdfs', 'MYC_neighborhood.rdf')

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
    with open(test_rdf_nfkb, 'rt') as fh:
        rdf_str_nfkb = fh.read()
    bp = bel.process_belrdf(rdf_str_nfkb)
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

#rdf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                        '..', '..', 'data')

def test_get_transcription():
    #with open(os.path.join(rdf_path, 'myc_neighborhood.rdf')) as f:
    #    myc_str = f.read()
    #    myc_bp = bel.process_belrdf(myc_str)
    #with open(test_rdf_myc, 'rt') as fh:
    #    rdf_str_myc = fh.read()
    #bp = bel.process_belrdf(rdf_str_myc)
    #transcription_stmts = []
    #for stmt in bp.statements + bp.indirect_stmts:
    #   if isinstance(stmt, RegulateAmount):
    #        transcription_stmts.append(stmt)
    #assert len(transcription_stmts) == 8
    pass

