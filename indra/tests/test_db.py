
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
    pmid = u'1234'
    text_ref_id = db.insert_text_ref(pmid=pmid)
    entries = db.get_text_refs_by_pmid(pmid)
    assert_equal(len(entries), 1, "One thing inserted, multiple entries found.")
    assert_equal(entries[0].pmid, pmid)#, "Got back the wrong pmid.")
    assert_equal(entries[0].id, text_ref_id, "Got back wrong text_ref_id.")
    db._clear()

def test_uniqueness_text_ref_doi_pmid():
    "Test uniqueness enforcement behavior for text_ref insertion."
    db = startup()
    pmid = u'1234'
    doi = u'foo/1234'
    text_ref_id_1 = db.insert_text_ref(doi=doi, pmid=pmid)
    text_ref_id_2 = db.insert_text_ref(doi=doi, pmid=pmid)
    entries = db.get_text_refs_by_pmid(pmid)
    assert_equal(len(entries), 1, 'Got multiple entries, no uniqueness.')
    assert_equal(text_ref_id_1, text_ref_id_2)
    db._clear()

def test_uniqueness_text_ref_url():
    "Test whether the uniqueness imposed on the url of text_refs is enforced."
    db = startup()
    url = u'http://foobar.com'
    db.insert_text_ref(url=url)
    try:
        db.insert_text_ref(url=url)
    except IntegrityError:
        return # PASS
    assert False, "Uniqueness was not enforced."
    
    
    