import logging
import unittest
from nose.plugins.attrib import attr
from indra.literature.adeft_tools import universal_extract_paragraphs, \
    filter_paragraphs
from indra.literature import pmc_client, elsevier_client, pubmed_client

logger = logging.getLogger(__name__)


@attr('nonpublic', 'webservice')
@unittest.skip('Elsevier credentials currently not operational')
def test_universal_extract_paragraphs_elsevier():
    doi = '10.1016/B978-0-12-416673-8.00004-6'
    xml_str = elsevier_client.download_article(doi)
    paragraphs = universal_extract_paragraphs(xml_str)
    if len(paragraphs) <= 1:
        logger.warning('Unable to extract paragraphs from XML string:\n'
                       '%s...' % xml_str[:2000])
    assert len(paragraphs) > 1


@attr('webservice')
def test_universal_extract_paragraphs_pmc():
    pmc_id = 'PMC3262597'
    xml_str = pmc_client.get_xml(pmc_id)
    paragraphs = universal_extract_paragraphs(xml_str)
    assert len(paragraphs) > 1


@attr('webservice')
def test_universal_extract_paragraphs_abstract():
    pmid = '16511588'
    abstract = pubmed_client.get_abstract(pmid)
    result = universal_extract_paragraphs(abstract)
    assert result[0] == abstract


def test_universal_extract_texts_contains():
    example = ['eeeeeeeeeeeeeEReeeeeeeeee',
               'eeeeeee-ER-eeeeeeeeeeeeee',
               'eeeeeee ER eeeeeeeeeeeeee',
               'eeeeeeeER eeeeeeeeeeeeeee']
    result = ('eeeeeee-ER-eeeeeeeeeeeeee\n'
              'eeeeeee ER eeeeeeeeeeeeee\n')
    text = filter_paragraphs(example, contains='ER')
    assert text == result


def test_universal_extract_texts_contains_union():
    example = ['eeeeeeeeeeNPeeeeeeeeeeeeee',
               'eeeeeee-NPs-eeeeeeeeeeeeee',
               'eeeeeee NP eeeeeeeeeeeeeee',
               'eeeeeeeNPseeeee NP-eeeeeee',
               'eeeeeeeeeeeeeeeeeeeeeeeeee']
    result = ('eeeeeee-NPs-eeeeeeeeeeeeee\n'
              'eeeeeee NP eeeeeeeeeeeeeee\n'
              'eeeeeeeNPseeeee NP-eeeeeee\n')
    text = filter_paragraphs(example, contains=['NP', 'NPs'])
    assert text == result
