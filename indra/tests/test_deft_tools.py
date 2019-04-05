from nose.plugins.attrib import attr

from indra.literature.deft_tools import universal_extract_text
from indra.literature import pmc_client, elsevier_client, pubmed_client


@attr('nonpublic', 'webservice')
def test_universal_extract_text_elsevier():
    doi = '10.1016/B978-0-12-416673-8.00004-6'
    xml_str = elsevier_client.download_article(doi)
    text = universal_extract_text(xml_str)
    assert text is not None
    assert ' ER ' in text


@attr('webservice')
def test_universal_extract_text_pmc():
    pmc_id = 'PMC3262597'
    xml_str = pmc_client.get_xml(pmc_id)
    text = universal_extract_text(xml_str)
    assert text is not None
    assert ' ER ' in text


@attr('webservice')
def test_universal_extract_text_abstract():
    pmid = '16511588'
    abstract = pubmed_client.get_abstract(pmid)
    result = universal_extract_text(abstract)
    assert result == abstract + '\n'


@attr('webservice')
def test_universal_extract_text_contains():
    pmc_id = 'PMC3262597'
    xml_str = pmc_client.get_xml(pmc_id)
    text1 = universal_extract_text(xml_str)
    text2 = universal_extract_text(xml_str, contains='ER')
    assert text1 is not None
    assert text2 is not None
    assert ' ER ' in text1 and ' ER ' in text2
    assert len(text2) < len(text1)


@attr('webservice')
def test_universal_extract_text_contains_union():
    pmc_id = 'PMC4954987'
    xml_str = pmc_client.get_xml(pmc_id)
    text1 = universal_extract_text(xml_str)
    text2 = universal_extract_text(xml_str, contains='NP')
    text3 = universal_extract_text(xml_str, contains='NPs')
    text4 = universal_extract_text(xml_str, contains=['NP', 'NPs'])
    assert text1 is not None
    assert text2 is not None
    assert text3 is not None
    assert text4 is not None
    assert len(text2) < len(text3) < len(text4) < len(text1)
