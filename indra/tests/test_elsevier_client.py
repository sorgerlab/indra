from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.literature import elsevier_client as ec

def test_get_fulltext_article():
    # This article is not open access so in order to get a full text response
    # with a body element requires full text access keys to be correctly
    # set up.
    doi = '10.1016/j.cell.2016.02.059'
    text = ec.get_article(doi)
    assert text is not None

def test_get_abstract():
    # If we have an API key but are not on an approved IP or don't have a
    # necessary institution key, we should still be able to get the abstract.
    # If there is a problem with the API key itself, this will log and error
    # and return None.
    doi = '10.1016/j.cell.2016.02.059'
    text = ec.get_abstract(doi)
    assert text is not None

def test_get_converted_article_body():
    """Make sure we can get fulltext of an article that has
    ja:converted-article as its principal sub-element."""
    # PMID: 11851341
    doi = '10.1006/jmbi.2001.5334'
    xml_str = ec.download_article(doi)
    body = ec.extract_text(xml_str)
    assert body

def test_get_rawtext():
    """Make sure we can get content of an article that has content in
    xocs:rawtext"""
    # PMID: 20072652
    doi = '10.1593/neo.91196'
    xml_str = ec.download_article(doi)
    body = ec.extract_text(xml_str)
    assert body

def test_article():
    # PMID: 11302724
    doi = '10.1006/bbrc.2001.4693'
    xml_str = ec.download_article(doi)
    body = ec.extract_text(xml_str)
    assert body is None
