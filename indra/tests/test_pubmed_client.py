import time
from indra.literature import pubmed_client
from nose.plugins.attrib import attr


@attr('webservice')
def test_get_ids():
    time.sleep(0.5)
    ids = pubmed_client.get_ids('braf', retmax=10, db='pubmed')
    assert len(ids) == 10


@attr('webservice')
def test_get_no_ids():
    time.sleep(0.5)
    ids = pubmed_client.get_ids('UUuXNWMCusRpcVTX', retmax=10, db='pubmed')
    assert not ids


@attr('webservice')
def test_get_ids():
    time.sleep(0.5)
    ids1 = pubmed_client.get_ids('JUN', use_text_word=False)
    ids2 = pubmed_client.get_ids('JUN', use_text_word=True)
    assert len(ids1) > len(ids2)


@attr('webservice')
def test_get_id_count():
    time.sleep(0.5)
    id_count = pubmed_client.get_id_count('SDLFKJSLDKJH')
    assert id_count == 0
    id_count = pubmed_client.get_id_count('KRAS')
    assert id_count > 0


@attr('webservice')
def test_get_id_mesh():
    time.sleep(0.5)
    ids = pubmed_client.get_ids_for_mesh('D009101')
    assert len(ids) > 35000
    ids_maj = pubmed_client.get_ids_for_mesh('D009101', major_topic=True)
    assert len(ids_maj) < len(ids)


@attr('webservice')
def test_get_pmc_ids():
    time.sleep(0.5)
    ids = pubmed_client.get_ids('braf', retmax=10, db='pmc')
    assert len(ids) == 10
    assert all(int(i[0]) >= 5 for i in ids), ids


@attr('webservice')
def test_get_title():
    time.sleep(0.5)
    title = pubmed_client.get_title('27754804')
    assert title
    assert title.lower().startswith('targeting autophagy')


@attr('webservice')
def test_get_title_prefix():
    time.sleep(0.5)
    title = pubmed_client.get_title('PMID27754804')
    assert title
    assert title.lower().startswith('targeting autophagy')


@attr('webservice')
def test_expand_pagination():
    time.sleep(0.5)
    pages = '456-7'
    new_pages = pubmed_client.expand_pagination(pages)
    assert new_pages == '456-457'


@attr('webservice')
def test_get_abstract_notitle():
    time.sleep(0.5)
    abstract = pubmed_client.get_abstract('27754804', prepend_title=False)
    assert abstract.startswith('The RAF inhibitor')
    assert abstract.endswith('vemurafenib.')


@attr('webservice')
def test_get_abstract_title():
    time.sleep(0.5)
    abstract = pubmed_client.get_abstract('27754804', prepend_title=True)
    assert abstract.lower().startswith('targeting autophagy')
    assert abstract.endswith('vemurafenib.')


@attr('webservice')
def test_get_abstract2():
    time.sleep(0.5)
    # Try another one
    abstract = pubmed_client.get_abstract('27123883')


@attr('webservice')
def test_get_no_abstract():
    time.sleep(0.5)
    abstract = pubmed_client.get_abstract('xx')
    assert abstract is None


@attr('webservice')
def test_get_ids_for_gene():
    time.sleep(0.5)
    ids = pubmed_client.get_ids_for_gene('EXOC1')
    assert ids


@attr('webservice')
def test_get_metadata_for_ids():
    time.sleep(0.5)
    pmids = ['27123883', '27121204', '27115606']
    metadata = pubmed_client.get_metadata_for_ids(pmids)


@attr('webservice')
def test_get_pub_date():
    time.sleep(0.5)
    pmids = ['27123883', '27121204', '27115606']
    metadata = pubmed_client.get_metadata_for_ids(pmids)
    assert metadata[pmids[0]]['publication_date']['year'] == 2016
    assert metadata[pmids[0]]['publication_date']['month'] == 4
    assert metadata[pmids[0]]['publication_date']['day'] == 29
    assert metadata[pmids[1]]['publication_date']['year'] == 2016
    assert metadata[pmids[1]]['publication_date']['month'] == 4
    assert metadata[pmids[1]]['publication_date']['day'] == 29
    assert metadata[pmids[2]]['publication_date']['year'] == 2016
    assert metadata[pmids[2]]['publication_date']['month'] == 4
    assert metadata[pmids[2]]['publication_date']['day'] == 27


@attr('webservice')
def test_send_request_invalid():
    time.sleep(0.5)
    res = pubmed_client.send_request('http://xxxxxxx', data={})
    assert res is None


@attr('webservice')
def test_abstract_with_html_embedded():
    time.sleep(0.5)
    res = pubmed_client.get_abstract('25484845')
    assert len(res) > 4, res


@attr('webservice')
def test_pmid_27821631():
    time.sleep(0.5)
    pmid = '27821631'
    res = pubmed_client.get_abstract(pmid)
    assert len(res) > 50, res
    res = pubmed_client.get_metadata_for_ids([pmid], get_abstracts=True)
    assert res[pmid]['title'] is not None
    assert len(res[pmid]['abstract']) > 50


@attr('webservice')
def test_get_annotations():
    time.sleep(0.5)
    pmid = '30971'
    tree = pubmed_client.send_request(pubmed_client.pubmed_fetch,
                                      dict(db='pubmed', retmode='xml',
                                           id=pmid))
    results = pubmed_client.get_metadata_from_xml_tree(tree,
                                                       mesh_annotations=True)
    assert len(results) == 1, len(results)
    assert 'mesh_annotations' in results[pmid].keys(), results[pmid].keys()
    me_ans = results[pmid]['mesh_annotations']
    assert len(me_ans) == 9, len(me_ans)
    assert all(d['mesh'].startswith('D') for d in me_ans)
    assert any(d['major_topic'] for d in me_ans)

