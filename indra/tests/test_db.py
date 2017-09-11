from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from os import listdir, remove
from nose import SkipTest
from nose.tools import assert_equal
from sqlalchemy.exc import IntegrityError
from indra.db.populate_content import Medline, PmcOA, Manuscripts
from indra.db import DatabaseManager, get_aws_db, texttypes

TEST_FILE = 'indra_test.db'
TEST_HOST = 'sqlite:///' + TEST_FILE

#TODO: implement setup-teardown system.
LOCAL_START_SUCCESS = True
REMOTE_START_SUCCESS = True

#==============================================================================
# The following are some helpful functions for the rest of the tests.
#==============================================================================
def assert_contents_equal(list1, list2, msg = None):
    "Check that the contenst of two lists are the same, regardless of order."
    res = set(list1) == set(list2)
    err_msg = "Contents of lists do not match:\n%s\n%s\n" % (list1, list2)
    if msg is not None:
        err_msg += msg
    assert res, err_msg


def local_startup():
    "Set up the database for testing."
    if not LOCAL_START_SUCCESS:
        raise SkipTest("Could not create test table.")
    
    db = DatabaseManager(TEST_HOST, sqltype='sqlite')
    db.grab_session()
    db._clear()
    
    return db

def remote_startup():
    "Set up the database for testing."
    if not REMOTE_START_SUCCESS:
        raise SkipTest("Could not create test table.")
    
    db = get_aws_db()
    db.grab_session()
    db._clear()
    
    return db

#==============================================================================
# The following are tests for the database manager itself.
#==============================================================================
def test_local_startup():
    "Test the local_startup function."
    # Cleanup from the last run, if necessary.
    if TEST_FILE in listdir('.'):
        remove(TEST_FILE)
    try:
        db = local_startup()
        db.create_tables()
        assert TEST_FILE in listdir('.'), "Test database not created"
        assert db.session is not None, "Could not get db session."
    except:
        global LOCAL_START_SUCCESS
        LOCAL_START_SUCCESS = False
        raise
    return


def test_remote_startup():
    "Test the remote_startup function."
    # Cleanup from the last run, if necessary.
    try:
        db = remote_startup()
        db.create_tables()
        assert db.session is not None, "Could not get db session."
    except:
        global REMOTE_START_SUCCESS
        REMOTE_START_SUCCESS = False
        raise
    return

def test_create_tables():
    "Test the create_tables feature"
    db = local_startup()
    db.create_tables()
    assert_contents_equal(db.get_active_tables(), db.get_tables())


def test_insert_and_query_pmid():
    "Test that we can add a text_ref and get the text_ref back."
    db = local_startup()
    pmid = '1234'
    text_ref_id = db.insert('text_ref', pmid=pmid)
    entries = db.select('text_ref', db.TextRef.pmid==pmid)
    assert_equal(len(entries), 1, "One thing inserted, multiple entries found.")
    assert_equal(entries[0].pmid, pmid)#, "Got back the wrong pmid.")
    assert_equal(entries[0].id, text_ref_id, "Got back wrong text_ref_id.")


def test_uniqueness_text_ref_doi_pmid():
    "Test uniqueness enforcement behavior for text_ref insertion."
    db = local_startup()
    pmid = '1234'
    doi = 'foo/1234'
    db.insert('text_ref', doi=doi, pmid=pmid)
    try:
        db.insert('text_ref', doi=doi, pmid=pmid)
    except IntegrityError:
        return # PASS
    finally:
        db._clear()
    assert False, "Uniqueness was not enforced."


def test_uniqueness_text_ref_url():
    "Test whether the uniqueness imposed on the url of text_refs is enforced."
    db = local_startup()
    url = 'http://foobar.com'
    db.insert('text_ref', url=url)
    try:
        db.insert('text_ref', url=url)
    except IntegrityError:
        return # PASS
    assert False, "Uniqueness was not enforced."


def test_get_abstracts():
    "Test the ability to get a list of abstracts."
    db = local_startup()
    
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
    not_found_fmt = 'If found, something is probably wrong with %s.'
    db.insert_many(
        'text_content',
        [
            {
                'text_ref_id':ref_id_list[0], #pmid=1234
                'source':'God',
                'format':'stone tablet',
                'text_type':'abstract',
                'content':(found_abst_fmt % '1234').encode('utf8')
                },
            {
                'text_ref_id':ref_id_list[1], #pmid=5678
                'source':'Satan',
                'format':'blood',
                'text_type':'full_content',
                'content':(not_found_fmt % 'text_type filter').encode('utf8')
                },
            {
                'text_ref_id':ref_id_list[1], #pmid=5678
                'source':'Satan',
                'format':'blood',
                'text_type':'abstract',
                'content':(found_abst_fmt % '5678').encode('utf8')
                },
            {
                'text_ref_id':ref_id_list[2], #no pmid
                'source':'Nature',
                'format':'tears',
                'text_type':'abstract',
                'content':(not_found_fmt % 'text_ref_id filter').encode('utf8')
                },
            {
                'text_ref_id':ref_id_list[3], #pmid=1357
                'source':'A Voice Inside Your Head',
                'format':'whispers',
                'text_type':'abstract',
                'content':(not_found_fmt % 'pmid filter').encode('utf8')
                }
            ]
        )
    
    expected = [(pmid, (found_abst_fmt % pmid).encode('utf8')) for pmid in ['1234', '5678']]
    received = db.get_abstracts_by_pmids(['1234', '5678', '1357'], unzip=False)
    assert_contents_equal(expected, received, "Did not get expected abstracts.")


def test_get_all_pmids():
    "Test whether we get all the pmids."
    db = local_startup()
    db.insert_many('text_ref', [{'pmid':'1234'}, {'pmid':'5678'}])
    pmid_list = db.get_all_pmids()
    assert_contents_equal(pmid_list, ['1234','5678'])


#==============================================================================
# The following are tests for the initial population of the database. This
# includes uploading data from Medline, PMC, Springer, and Elsevier. These
# tend to make greater use of the database, and are likely to be slower.
#==============================================================================
def test_full_local_upload():
    "Test whether we can perform a targeted upload to a local db."
    # This uses a specially curated sample directory designed to access all code
    # paths that the real system might experience, but on a much smaller (thus
    # faster) scale. Errors in the ftp service will not be caught by this test.
    db = local_startup()
    loc_path = 'test_ftp'
    Medline(ftp_url=loc_path, local=True).populate(db)
    tr_list = db.select('text_ref')
    assert len(tr_list), "No text refs were added..."
    assert all([hasattr(tr, 'pmid') for tr in tr_list]),\
        'All text_refs MUST have pmids by now.'
    #assert all([hasattr(tr, 'pmcid') for tr in tr_list]),\
    #    'All text_refs should have pmcids at this point.' 
    PmcOA(ftp_url=loc_path, local=True).populate(db)
    tc_list = db.select('text_content', db.TextContent.text_type==texttypes.FULLTEXT)
    assert len(tc_list), "No fulltext was added."
    Manuscripts(ftp_url=loc_path, local=True).populate(db)
    tc_list = db.select('text_content', db.TextContent.source==Manuscripts.my_source)
    assert len(tc_list), "No manuscripts uploaded."


def test_full_remote_upload():
    "Test whether we can perform a targeted upload to aws."
    # This uses a specially curated sample directory designed to access all code
    # paths that the real system might experience, but on a much smaller (thus
    # faster) scale. Errors in the ftp service will not be caught by this test.
    db = remote_startup()
    loc_path = 'test_ftp'
    Medline(ftp_url=loc_path, local=True).populate(db)
    tr_list = db.select('text_ref')
    assert all([hasattr(tr, 'pmid') for tr in tr_list]),\
        'All text_refs MUST have pmids by now.'
    #assert all([hasattr(tr, 'pmcid') for tr in tr_list]),\
    #    'All text_refs should have pmcids at this point.' 
    PmcOA(ftp_url=loc_path, local=True).populate(db)
    tc_list = db.select('text_content', db.TextContent.text_type==texttypes.FULLTEXT)
    assert len(tc_list), "No fulltext was added."
    Manuscripts(ftp_url=loc_path, local=True).populate(db)
    tc_list = db.select('text_content', db.TextContent.source==Manuscripts.my_source)
    assert len(tc_list), "No manuscripts uploaded."


def test_ftp_service():
    "Test the progenitor childrens' ftp access."
    cases = [
        ('.csv', 'csv_as_dict'), 
        #('.xml.gz', 'xml_file'),
        ('.txt', 'file')
        ]
    for Child in [Medline, PmcOA, Manuscripts]:
        c = Child()
        files = c.ftp_ls()
        for end, meth in cases:
            rel_files = [f for f in files if f.endswith(end)]
            # TODO: check metadata to choose small files.
            if len(rel_files) > 0:
                file = rel_files[0]
                ret = getattr(c, 'get_' + meth)(file)
                assert ret is not None,\
                    "Failed to load %s from %s." % (file, Child.__name__)
        # TODO: test the download.



