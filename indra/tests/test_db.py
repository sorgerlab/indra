from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from os import listdir, remove
from indra.db import DatabaseManager
from nose import SkipTest
from nose.tools import assert_equal
from sqlalchemy.exc import IntegrityError

TEST_FILE = 'indra_test.db'
TEST_HOST = 'sqlite:///' + TEST_FILE

#TODO: implement setup-teardown system.
START_SUCCESS = True

def assert_contents_equal(list1, list2, msg = None):
    "Check that the contenst of two lists are the same, regardless of order."
    res = set(list1) == set(list2)
    err_msg = "Contents of lists do not match:\n%s\n%s\n" % (list1, list2)
    if msg is not None:
        err_msg += msg
    assert res, err_msg


def startup():
    "Set up the database for testing."
    if not START_SUCCESS:
        raise SkipTest("Could not create test table.")
    
    db = DatabaseManager(TEST_HOST, sqltype='sqlite')
    db.get_session()
    db._clear()
    
    return db


def test_startup():
    "Test the startup function."
    # Cleanup from the last run, if necessary.
    if TEST_FILE in listdir('.'):
        remove(TEST_FILE)
    try:
        db = startup()
        db.create_tables()
        assert TEST_FILE in listdir('.'), "Test database not created"
        assert db.session is not None, "Could not get db session."
    except:
        global START_SUCCESS
        START_SUCCESS = False
        raise
    return

def test_create_tables():
    "Test the create_tables feature"
    db = startup()
    db.create_tables()
    assert_contents_equal(db.get_active_tables(), db.get_tables())


def test_insert_and_query_pmid():
    "Test that we can add a text_ref and get the text_ref back."
    db = startup()
    pmid = '1234'
    text_ref_id = db.insert_text_ref(pmid=pmid)
    entries = db.get_text_refs_by_pmid(pmid)
    assert_equal(len(entries), 1, "One thing inserted, multiple entries found.")
    assert_equal(entries[0].pmid, pmid)#, "Got back the wrong pmid.")
    assert_equal(entries[0].id, text_ref_id, "Got back wrong text_ref_id.")


def test_uniqueness_text_ref_doi_pmid():
    "Test uniqueness enforcement behavior for text_ref insertion."
    db = startup()
    pmid = '1234'
    doi = 'foo/1234'
    db.insert_text_ref(doi=doi, pmid=pmid)
    try:
        db.insert_text_ref(doi=doi, pmid=pmid)
    except IntegrityError:
        return # PASS
    finally:
        db._clear()
    assert False, "Uniqueness was not enforced."


def test_uniqueness_text_ref_url():
    "Test whether the uniqueness imposed on the url of text_refs is enforced."
    db = startup()
    url = 'http://foobar.com'
    db.insert_text_ref(url=url)
    try:
        db.insert_text_ref(url=url)
    except IntegrityError:
        return # PASS
    assert False, "Uniqueness was not enforced."


def test_get_abstracts():
    "Test the ability to get a list of abstracts."
    db = startup()
    
    # Create a world of abstracts.
    ref_id_list = db.insert_many(
        'text_ref',
        [
            {'pmid':'1234'}, # searched for, just abstract.
            {'pmid':'5678'}, # searched for, abstract and full_text
            {'doi':'foo/234'}, # content, but no pmid.
            {'pmid':'2468'}, # not searched for.
            {'pmid':'1357'} # searched for, but no conent.
            ]
        )
    found_abst_fmt = 'This should be found alongside pmid %s.'
    not_found_fmt = 'If found, something is wrong with %s.'
    db.insert_many(
        'text_content',
        [
            {
                'text_ref_id':ref_id_list[0], #pmid=1234
                'source':'God',
                'format':'stone tablet',
                'text_type':'abstract',
                'content':found_abst_fmt % '1234'
                },
            {
                'text_ref_id':ref_id_list[1], #pmid=5678
                'source':'Satan',
                'format':'blood',
                'text_type':'full_content',
                'content':not_found_fmt % 'text_type filter'
                },
            {
                'text_ref_id':ref_id_list[1], #pmid=5678
                'source':'Satan',
                'format':'blood',
                'text_type':'abstract',
                'content':found_abst_fmt % '5678'
                },
            {
                'text_ref_id':ref_id_list[2], #no pmid
                'source':'Nature',
                'format':'tears',
                'text_type':'abstract',
                'content':not_found_fmt % 'text_ref_id filter'
                },
            {
                'text_ref_id':ref_id_list[3], #pmid=1357
                'source':'A Voice Inside Your Head',
                'format':'whispers',
                'text_type':'abstract',
                'content':not_found_fmt % 'pmid filter.'
                }
            ]
        )
    
    expected = [(pmid, found_abst_fmt % pmid) for pmid in ['1234', '5678']]
    received = db.get_abstracts_by_pmids(['1234', '5678', '1357'], unzip=False)
    assert_contents_equal(expected, received, "Did not get expected abstracts.")

def test_get_all_pmids():
    "Test whether we get all the pmids."
    db = startup()
    db.insert_many('text_ref', [{'pmid':'1234'}, {'pmid':'5678'}])
    pmid_list = db.get_all_pmids()
    assert_contents_equal(pmid_list, ['1234','5678'])

