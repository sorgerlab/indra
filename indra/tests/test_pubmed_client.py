from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.literature import pubmed_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr


@attr('webservice')
def test_get_ids():
    ids = pubmed_client.get_ids('braf', retmax=10, db='pubmed')
    assert(len(ids) == 10)
    assert unicode_strs(ids)


@attr('webservice')
def test_get_no_ids():
    ids = pubmed_client.get_ids('xkcd', retmax=10, db='pubmed')
    assert(not ids)


@attr('webservice')
def test_get_ids():
    ids1 = pubmed_client.get_ids('JUN', use_text_word=False)
    ids2 = pubmed_client.get_ids('JUN', use_text_word=True)
    assert(len(ids1) > len(ids2))
    assert unicode_strs(ids1)
    assert unicode_strs(ids2)


@attr('webservice')
def test_get_pmc_ids():
    ids = pubmed_client.get_ids('braf', retmax=10, db='pmc')
    assert(len(ids) == 10)
    assert(len([i for i in ids if i.startswith('5') or
                i.startswith('4')]) == 10)
    assert unicode_strs(ids)


@attr('webservice')
def test_get_title():
    title = pubmed_client.get_title('27754804')
    assert(title)
    assert(title.lower().startswith('targeting autophagy'))


@attr('webservice')
def test_get_title_prefix():
    title = pubmed_client.get_title('PMID27754804')
    assert(title)
    assert(title.lower().startswith('targeting autophagy'))


@attr('webservice')
def test_expand_pagination():
    pages = '456-7'
    new_pages = pubmed_client.expand_pagination(pages)
    assert(new_pages == '456-457')


@attr('webservice')
def test_get_abstract_notitle():
    abstract = pubmed_client.get_abstract('27754804', prepend_title=False)
    assert(abstract.startswith('The RAF inhibitor'))
    assert(abstract.endswith('vemurafenib.'))
    assert unicode_strs(abstract)


@attr('webservice')
def test_get_abstract_title():
    abstract = pubmed_client.get_abstract('27754804', prepend_title=True)
    assert(abstract.lower().startswith('targeting autophagy'))
    assert(abstract.endswith('vemurafenib.'))
    assert unicode_strs(abstract)


@attr('webservice')
def test_get_abstract2():
    # Try another one
    abstract = pubmed_client.get_abstract('27123883')
    assert unicode_strs(abstract)


@attr('webservice')
def test_get_no_abstract():
    abstract = pubmed_client.get_abstract('xx')
    assert(abstract is None)


@attr('webservice')
def test_get_ids_for_gene():
    ids = pubmed_client.get_ids_for_gene('EXOC1')
    assert ids
    assert unicode_strs(ids)


@attr('webservice')
def test_get_metadata_for_ids():
    pmids = ['27123883', '27121204', '27115606']
    metadata = pubmed_client.get_metadata_for_ids(pmids)
    assert unicode_strs(metadata)


@attr('webservice')
def test_send_request_invalid():
    res = pubmed_client.send_request('http://xxxxxxx', data={})
    assert res is None


@attr('webservice')
def test_abstract_with_html_embedded():
    res = pubmed_client.get_abstract('25484845')
    assert len(res) > 4, res
