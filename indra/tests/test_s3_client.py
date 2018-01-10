from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.literature import s3_client
from indra.util import unicode_strs
import zlib
from nose.plugins.attrib import attr


@attr('webservice', 'nonpublic')
def test_check_pmid():
    pmid = s3_client.check_pmid(12345)
    assert pmid == 'PMID12345'
    assert unicode_strs(pmid)
    pmid = s3_client.check_pmid('12345')
    assert pmid == 'PMID12345'
    assert unicode_strs(pmid)
    pmid = s3_client.check_pmid('PMID12345')
    assert pmid == 'PMID12345'
    assert unicode_strs(pmid)


@attr('webservice', 'nonpublic')
def test_get_pmid_key():
    pmid = '12345'
    pmid_key = s3_client.get_pmid_key(pmid)
    assert pmid_key == s3_client.prefix + 'PMID12345'
    assert unicode_strs(pmid_key)


@attr('webservice', 'nonpublic')
def test_filter_keys():
    pmid_key = s3_client.get_pmid_key('1001287')
    key_list = s3_client.filter_keys(pmid_key)
    assert len(key_list) == 4


@attr('webservice', 'nonpublic')
def test_get_gz_object():
    # Get XML
    key = 'papers/PMID27297883/fulltext/txt'
    obj = s3_client.get_gz_object(key)
    assert unicode_strs(obj)
    # Get reach output
    key = 'papers/PMID27297883/reach'
    obj = s3_client.get_gz_object(key)
    assert unicode_strs(obj)


@attr('webservice', 'nonpublic')
def test_get_gz_object_nosuchkey():
    obj = s3_client.get_gz_object('foobar')
    assert obj is None


@attr('webservice', 'nonpublic')
def test_get_full_text():
    (content, content_type) = s3_client.get_full_text('27297883')
    assert unicode_strs((content, content_type))
    assert content_type == 'txt'
    (content, content_type) = s3_client.get_full_text('1001287')
    assert unicode_strs((content, content_type))
    assert content_type == 'pmc_oa_xml'
    # TODO: Find a paper that has only abstract
    #(content, content_type) = s3_client.get_full_text('27653174')
    #assert unicode_strs((content, content_type))
    #assert content_type == 'abstract'
    (content, content_type) = s3_client.get_full_text('000000')
    assert content is None and content_type is None


@attr('webservice', 'nonpublic')
def test_put_full_text():
    full_text = 'test_put_full_text'
    pmid_test = 'PMID000test1'
    s3_client.put_full_text(pmid_test, full_text, full_text_type='pmc_oa_txt')
    # Get the full text back
    (content, content_type) = s3_client.get_full_text(pmid_test)
    assert content == full_text
    assert content_type == 'pmc_oa_txt'
    assert unicode_strs(content)


@attr('webservice', 'nonpublic')
def test_put_abstract():
    abstract = 'test_put_abstract'
    pmid_test = 'PMID000test2'
    s3_client.put_abstract(pmid_test, abstract)
    # Get the abstract back
    (content, content_type) = s3_client.get_full_text(pmid_test)
    assert content == abstract
    assert content_type == 'abstract'
    assert unicode_strs(content)


@attr('webservice', 'nonpublic')
def test_reach_output():
    # Test put_reach_output
    reach_data = {'foo': 1, 'bar': {'baz': 2}}
    pmid = 'PMID000test3'
    reach_version = '42'
    source_text = 'pmc_oa_txt'
    s3_client.put_reader_output('reach', reach_data, pmid, reach_version, source_text)
    # Now get the data back
    retrieved_reach_data = s3_client.get_reader_output('reach', pmid)
    assert retrieved_reach_data == reach_data
    assert unicode_strs(retrieved_reach_data)
    # Get the reach version of the key we created
    ret_reach_version, ret_source_text = \
        s3_client.get_reader_metadata('reach', pmid)
    assert ret_reach_version == reach_version
    assert ret_source_text == source_text
    assert unicode_strs(ret_reach_version)


def test_gzip_string():
    content = 'asdf'
    content_enc = s3_client.gzip_string(content, 'content')
    content_dec = zlib.decompress(content_enc, 16+zlib.MAX_WBITS)
    content_dec_uni = content_dec.decode('utf-8')
    assert content == content_dec_uni


@attr('webservice', 'nonpublic')
def test_get_upload_content():
    pmid_s3_no_content = 'PMID000foobar'
    (ct, ct_type) = s3_client.get_upload_content(pmid_s3_no_content)
    assert ct is None
    assert ct_type is None

    pmid_s3_abstract_only = 'PMID000test4'
    s3_client.put_abstract(pmid_s3_abstract_only, 'foo')
    (ct, ct_type) = s3_client.get_upload_content(pmid_s3_abstract_only)
    assert ct == 'foo'
    assert ct_type == 'abstract'

    pmid_s3_fulltext = 'PMID000test5'
    s3_client.put_full_text(pmid_s3_fulltext, 'foo', full_text_type='txt')
    (ct, ct_type) = s3_client.get_upload_content(pmid_s3_fulltext)
    assert ct == 'foo'
    assert ct_type == 'txt'
