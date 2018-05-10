from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from rdflib.term import URIRef
from indra.util import unicode_strs
from indra.sources import bel
from indra.sources.bel.belrdf_processor import BelRdfProcessor
from indra.statements import RegulateAmount
from nose.plugins.attrib import attr

concept_prefix = 'http://www.openbel.org/bel/namespace//'
entity_prefix = 'http://www.openbel.org/bel/'

path_this = os.path.dirname(os.path.abspath(__file__))
test_rdf_nfkb = os.path.join(path_this, 'bel_rdfs', 'NFKB1_neighborhood.rdf')
test_rdf_myc = os.path.join(path_this, 'bel_rdfs', 'MYC_neighborhood.rdf')

def assert_pmids(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid is not None:
                assert ev.pmid.isdigit(), ev.pmid 

@attr('webservice')
def test_bel_ndex_query():
    bp = bel.process_ndex_neighborhood(['NFKB1'])
    assert_pmids(bp.statements)
    unicode_strs(bp.statements)


@attr('webservice')
def test_pybel_neighborhood_query():
    corpus = path_this + '/../../data/small_corpus.bel'
    bp = bel.process_pybel_neighborhood(['TP63'], corpus)
    assert bp.statements
    assert_pmids(bp.statements)
    unicode_strs(bp.statements)


def test_process_belrdf():
    with open(test_rdf_nfkb, 'rb') as fh:
        rdf_str_nfkb = fh.read().decode('utf-8')
    bp = bel.process_belrdf(rdf_str_nfkb)
    assert_pmids(bp.statements)
    unicode_strs(bp.statements)

def test_get_agent_up_from_hgnc():
    hgnc_sym = 'MAPK1'
    concept = concept_prefix + hgnc_sym
    entity = entity_prefix + 'p_HGNC_' + hgnc_sym
    ag = BelRdfProcessor._get_agent(concept, entity)
    assert ag.name == 'MAPK1'
    assert ag.db_refs.get('HGNC') == '6871'
    assert ag.db_refs.get('UP') == 'P28482'
    assert unicode_strs((concept, entity, ag))

def test_get_agent_hgnc_up_from_egid():
    entrez_id = '5594'
    concept = concept_prefix + entrez_id
    entity = entity_prefix + 'p_EGID_' + entrez_id
    ag = BelRdfProcessor._get_agent(concept, entity)
    assert ag.name == 'MAPK1'
    assert ag.db_refs.get('EGID') == entrez_id
    assert ag.db_refs.get('HGNC') == '6871'
    assert ag.db_refs.get('UP') == 'P28482'
    assert unicode_strs((concept, entity, ag))


def test_get_transcription():
    print("Opening file")
    with open(test_rdf_myc, 'rt') as fh:
        rdf_str_myc = fh.read()
    print("Process BEL RDF")
    bp = bel.process_belrdf(rdf_str_myc)
    transcription_stmts = []
    for stmt in bp.statements + bp.indirect_stmts:
       if isinstance(stmt, RegulateAmount):
            transcription_stmts.append(stmt)
    assert len(transcription_stmts) == 8
    pass

