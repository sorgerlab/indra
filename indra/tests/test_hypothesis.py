from nose.plugins.attrib import attr
from indra.sources import hypothesis
from indra.sources.hypothesis.processor import HypothesisProcessor, \
    get_context_entry, parse_grounding_entry, get_text_refs


@attr('nonpublic')
def test_process_indra_annnotations():
    hp = hypothesis.process_annotations()
    assert hp.statements
    for stmt in hp.statements:
        print(stmt)
        print(stmt.evidence[0])


def test_get_text_refs_pmid():
    url = 'https://www.ncbi.nlm.nih.gov/pubmed/32196952'
    refs = get_text_refs(url)
    assert refs.get('PMID') == '32196952', refs
    assert refs.get('URL') == url, refs


def test_get_text_refs_pmcid():
    url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7071777/'
    refs = get_text_refs(url)
    assert refs.get('PMCID') == 'PMC7071777', refs
    assert refs.get('URL') == url, refs
