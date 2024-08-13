import time
from indra.literature import pubmed_client
import pytest


@pytest.mark.webservice
def test_get_ids1():
    time.sleep(0.5)
    ids = pubmed_client.get_ids('braf', retmax=10, db='pubmed')
    assert len(ids) == 10


@pytest.mark.webservice
def test_get_no_ids():
    time.sleep(0.5)
    ids = pubmed_client.get_ids('UUuXNWMCusRpcVTX', retmax=10, db='pubmed')
    assert not ids


@pytest.mark.webservice
def test_get_ids2():
    time.sleep(0.5)
    ids1 = pubmed_client.get_ids('JUN', use_text_word=False, reldate=365)
    ids2 = pubmed_client.get_ids('JUN', use_text_word=True, reldate=365)
    assert len(ids1) > len(ids2)


@pytest.mark.webservice
def test_get_id_count():
    time.sleep(0.5)
    id_count = pubmed_client.get_id_count('SDLFKJSLDKJH')
    assert id_count == 0
    id_count = pubmed_client.get_id_count('KRAS')
    assert id_count > 0


@pytest.mark.webservice
def test_get_id_mesh():
    time.sleep(0.5)
    ids = pubmed_client.get_ids_for_mesh('D009101', reldate=365)
    assert len(ids) > 100, len(ids)
    ids_maj = pubmed_client.get_ids_for_mesh('D009101', major_topic=True,
                                             reldate=365)
    assert len(ids_maj) < len(ids)


@pytest.mark.webservice
def test_get_id_mesh_supc():
    time.sleep(0.5)
    ids = pubmed_client.get_ids_for_mesh('C000111')
    assert len(ids) > 8, len(ids)


@pytest.mark.webservice
def test_get_pmc_ids():
    time.sleep(0.5)
    ids = pubmed_client.get_ids('braf', retmax=10, db='pmc')
    assert len(ids) == 10


@pytest.mark.webservice
def test_get_title():
    time.sleep(0.5)
    title = pubmed_client.get_title('27754804')
    assert title
    assert title.lower().startswith('targeting autophagy')


@pytest.mark.webservice
def test_get_title_prefix():
    time.sleep(0.5)
    title = pubmed_client.get_title('PMID27754804')
    assert title
    assert title.lower().startswith('targeting autophagy')


@pytest.mark.webservice
def test_get_complex_title():
    time.sleep(0.5)
    title = pubmed_client.get_title('33463523')
    assert title
    assert title.lower().startswith('atomic structures')
    assert title.lower().endswith('vascular plants.')


@pytest.mark.webservice
def test_expand_pagination():
    time.sleep(0.5)
    pages = '456-7'
    new_pages = pubmed_client.expand_pagination(pages)
    assert new_pages == '456-457'


@pytest.mark.webservice
def test_get_abstract_notitle():
    time.sleep(0.5)
    abstract = pubmed_client.get_abstract('27754804', prepend_title=False)
    assert abstract.startswith('The RAF inhibitor')
    assert abstract.endswith('vemurafenib.')


@pytest.mark.webservice
def test_get_abstract_title():
    time.sleep(0.5)
    abstract = pubmed_client.get_abstract('27754804', prepend_title=True)
    assert abstract.lower().startswith('targeting autophagy')
    assert abstract.endswith('vemurafenib.')


@pytest.mark.webservice
def test_get_abstract2():
    time.sleep(0.5)
    # Try another one
    abstract = pubmed_client.get_abstract('27123883')


@pytest.mark.webservice
def test_get_no_abstract():
    time.sleep(0.5)
    abstract = pubmed_client.get_abstract('xx')
    assert abstract is None


@pytest.mark.webservice
def test_get_ids_for_gene():
    time.sleep(0.5)
    ids = pubmed_client.get_ids_for_gene('EXOC1')
    assert ids


@pytest.mark.webservice
def test_get_metadata_for_ids():
    time.sleep(0.5)
    pmids1 = ['27123883', '27121204', '27115606']
    metadata1 = pubmed_client.get_metadata_for_ids(pmids1)
    pmids2 = ['27123883']
    time.sleep(0.5)
    metadata2 = pubmed_client.get_metadata_for_ids(pmids2,
                                                   detailed_authors=True)
    assert all(isinstance(a, str) for a in metadata1[pmids1[0]]['authors'])
    assert all(isinstance(a, dict) for a in metadata2[pmids2[0]]['authors'])
    assert metadata1[pmids1[0]]['authors'][0] == 'Le Rhun'
    assert metadata2[pmids1[0]]['authors'][0]['last_name'] == 'Le Rhun'
    assert 'Lille' in \
        metadata2[pmids1[0]]['authors'][0]['affiliations'][0]['name']


@pytest.mark.webservice
def test_get_paper_references():
    time.sleep(0.5)
    pmids = ['27123883', '27121204', '27115606']
    test_pmid = '27121204'
    referenced_pmid = '25439075'
    metadata_1 = pubmed_client.get_metadata_for_ids(pmids, references_included='pmid')
    assert len(metadata_1[test_pmid]['references']) != 0
    assert metadata_1[test_pmid]['references'][0] == referenced_pmid

    metadata_2 = pubmed_client.get_metadata_for_ids(pmids, references_included='detailed')
    assert len(metadata_2[test_pmid]['references']) != 0
    assert metadata_2[test_pmid]['references'][0]['pmid'] == referenced_pmid


@pytest.mark.webservice
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


@pytest.mark.webservice
def test_send_request_invalid():
    time.sleep(0.5)
    res = pubmed_client.send_request('http://xxxxxxx', data={})
    assert res is None


@pytest.mark.webservice
def test_abstract_with_html_embedded():
    time.sleep(0.5)
    res = pubmed_client.get_abstract('25484845')
    assert len(res) > 4, res


@pytest.mark.webservice
def test_pmid_27821631():
    time.sleep(0.5)
    pmid = '27821631'
    res = pubmed_client.get_abstract(pmid)
    assert len(res) > 50, res
    res = pubmed_client.get_metadata_for_ids([pmid], get_abstracts=True)
    assert res[pmid]['title'] is not None
    assert len(res[pmid]['abstract']) > 50


@pytest.mark.webservice
def test_get_annotations():
    time.sleep(0.5)
    pmid = '30971'
    tree = pubmed_client.get_full_xml(pmid)
    results = pubmed_client.get_metadata_from_xml_tree(tree,
                                                       mesh_annotations=True)
    assert len(results) == 1, len(results)
    assert 'mesh_annotations' in results[pmid], results[pmid]
    me_ans = results[pmid]['mesh_annotations']
    assert len(me_ans) == 9, len(me_ans)
    assert all(d['mesh'].startswith('D') for d in me_ans)
    assert any(d['major_topic'] for d in me_ans)


@pytest.mark.webservice
def test_get_supplementary_annotations():
    time.sleep(0.5)
    pmid = '30105248'
    anns = pubmed_client.get_mesh_annotations(pmid)
    assert len(anns) == 7, anns
    assert anns[0]['type'] == 'main'
    assert anns[0]['mesh'] == 'D053839'
    assert len(anns[0]['qualifiers']) == 1
    assert anns[0]['qualifiers'][0] == anns[0]['qualifier']
    supp_ann = anns[-1]
    assert supp_ann['type'] == 'supplementary'
    assert supp_ann['mesh'] == 'C000623891'
    assert supp_ann['text'] == 'Tomato yellow leaf curl virus'


@pytest.mark.webservice
def test_get_substance_annotations():
    pmid = '27959613'
    mesh_ids = pubmed_client.get_substance_annotations(pmid)
    example_mesh_id = 'D009570'
    wrong_mesh_id = 'D0074447'
    assert example_mesh_id in mesh_ids
    assert wrong_mesh_id not in mesh_ids


def test_is_retracted():
    assert pubmed_client.is_retracted('35463694')
    assert not pubmed_client.is_retracted('36938926')
