from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from sys import version_info
from nose import SkipTest
from os import listdir, remove, path
from nose.tools import assert_equal, assert_list_equal
from sqlalchemy.exc import IntegrityError
from indra.db.populate_content import Medline, PmcOA, Manuscripts
from indra.db import DatabaseManager, get_primary_db, texttypes, DEFAULTS_FILE


TEST_FILE = 'indra_test.db'
TEST_HOST = 'sqlite:///' + TEST_FILE
with open(DEFAULTS_FILE, 'r') as f:
    pass  # TODO: Do this


# TODO: implement setup-teardown system.
LOCAL_START_SUCCESS = True
REMOTE_START_SUCCESS = False

#==============================================================================
# The following are some helpful functions for the rest of the tests.
#==============================================================================
def assert_contents_equal(list1, list2, msg=None):
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

    db = get_primary_db()
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
    except Exception:
        global LOCAL_START_SUCCESS
        LOCAL_START_SUCCESS = False
        raise
    return


def test_remote_startup():
    "Test the remote_startup function."
    # Cleanup from the last run, if necessary.
    raise SkipTest("Do not mess with the database...")
    try:
        db = remote_startup()
        db.create_tables()
        assert db.session is not None, "Could not get db session."
    except Exception:
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
    entries = db.select_all('text_ref', db.TextRef.pmid == pmid)
    assert_equal(len(entries), 1, "One item inserted, multiple entries found.")
    assert_equal(entries[0].pmid, pmid)
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
        return  # PASS
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
        return  # PASS
    assert False, "Uniqueness was not enforced."


def test_get_abstracts():
    "Test the ability to get a list of abstracts."
    db = local_startup()

    # Create a world of abstracts.
    ref_id_list = db.insert_many(
        'text_ref',
        [
            {'pmid': '1234'},  # searched for, just abstract.
            {'pmid': '5678'},  # searched for, abstract and full_text
            {'doi': 'foo/234'},  # content, but no pmid.
            {'pmid': '2468'},  # not searched for.
            {'pmid': '1357'}  # searched for, but no conent.
            ]
        )
    found_abst_fmt = 'This should be found alongside pmid %s.'
    not_found_fmt = 'If found, something is probably wrong with %s.'
    db.insert_many(
        'text_content',
        [
            {
                'text_ref_id': ref_id_list[0],  # pmid=1234
                'source':'God',
                'format':'stone tablet',
                'text_type':'abstract',
                'content':(found_abst_fmt % '1234').encode('utf8')
                },
            {
                'text_ref_id': ref_id_list[1],  # pmid=5678
                'source':'Satan',
                'format':'blood',
                'text_type':'full_content',
                'content':(not_found_fmt % 'text_type filter').encode('utf8')
                },
            {
                'text_ref_id': ref_id_list[1],  # pmid=5678
                'source':'Satan',
                'format':'blood',
                'text_type':'abstract',
                'content':(found_abst_fmt % '5678').encode('utf8')
                },
            {
                'text_ref_id': ref_id_list[2],  # no pmid
                'source':'Nature',
                'format':'tears',
                'text_type':'abstract',
                'content':(not_found_fmt % 'text_ref_id filter').encode('utf8')
                },
            {
                'text_ref_id': ref_id_list[3],  # pmid=1357
                'source':'A Voice Inside Your Head',
                'format':'whispers',
                'text_type':'abstract',
                'content':(not_found_fmt % 'pmid filter').encode('utf8')
                }
            ]
        )

    expected = [(pmid, (found_abst_fmt % pmid).encode('utf8'))
                for pmid in ['1234', '5678']]
    received = db.get_abstracts_by_pmids(['1234', '5678', '1357'], unzip=False)
    assert_contents_equal(expected, received, "Didn't get expected abstracts.")


def test_get_all_pmids():
    "Test whether we get all the pmids."
    db = local_startup()
    db.insert_many('text_ref', [{'pmid': '1234'}, {'pmid': '5678'}])
    pmid_list = db.get_all_pmids()
    assert_contents_equal(pmid_list, ['1234', '5678'])


#==============================================================================
# The following are tests for the initial population of the database. This
# includes uploading data from Medline, PMC, Springer, and Elsevier. These
# tend to make greater use of the database, and are likely to be slower.
#==============================================================================
TEST_POPULATE = True
if version_info.major is not 3:
    TEST_POPULATE = False

TEST_FTP = 'test_ftp'
if TEST_POPULATE and not path.exists(TEST_FTP):
    print("Creating test directory. This could take a while...")
    from indra.db.build_sample_set import build_set
    build_set(2, TEST_FTP)


def test_full_local_upload():
    "Test whether we can perform a targeted upload to a local db."
    # This uses a specially curated sample directory designed to access all code
    # paths that the real system might experience, but on a much smaller (thus
    # faster) scale. Errors in the ftp service will not be caught by this test.
    if not TEST_POPULATE:
        raise SkipTest("These feautures are only supported in Python 3.x")
    db = local_startup()
    Medline(ftp_url=TEST_FTP, local=True).populate(db)
    tr_list = db.select_all('text_ref')
    assert len(tr_list), "No text refs were added..."
    assert all([hasattr(tr, 'pmid') for tr in tr_list]),\
        'All text_refs MUST have pmids by now.'
    PmcOA(ftp_url=TEST_FTP, local=True).populate(db)
    tc_list = db.select_all(
        'text_content',
        db.TextContent.text_type == texttypes.FULLTEXT)
    assert len(tc_list), "No fulltext was added."
    Manuscripts(ftp_url=TEST_FTP, local=True).populate(db)
    tc_list = db.select_all(
        'text_content',
        db.TextContent.source == Manuscripts.my_source
        )
    assert len(tc_list), "No manuscripts uploaded."


def test_multiple_pmids():
    "Test that pre-existing pmids are correctly handled."
    if not TEST_POPULATE:
        raise SkipTest("These feautures are only supported in Python 3.x")
    db = local_startup()
    med = Medline(ftp_url=TEST_FTP, local=True)
    med.populate(db)
    num_refs = len(db.select_all('text_ref'))
    med.populate(db)
    assert len(db.select_all('text_ref')) == num_refs,\
        "Duplicate pmids allowed to be submitted.."
    return


def test_multible_pmc_oa_content():
    "Test to make sure repeated content is handled correctly."
    if not TEST_POPULATE:
        raise SkipTest("These feautures are only supported in Python 3.x")
    db = local_startup()
    pmc = PmcOA(ftp_url=TEST_FTP, local=True)
    pmc.populate(db)
    num_conts = len(db.select_all('text_content'))
    pmc.populate(db)
    assert len(db.select_all('text_content')) == num_conts,\
        "Duplicate text content allowed to be submitted."
    return


def test_multiple_text_ref_pmc_oa():
    "Test whether a duplicate text ref in pmc oa is handled correctly."
    if not TEST_POPULATE:
        raise SkipTest("These feautures are only supported in Python 3.x")
    db = local_startup()
    pmc = PmcOA(ftp_url=TEST_FTP, local=True)
    inp = dict.fromkeys(pmc.tr_cols)
    inp.update(pmcid='PMC5579538', doi='10.1021/acsomega.7b00205')
    pmc.upload_batch(db, [inp], [])
    num_refs = len(db.select_all('text_ref'))
    pmc.upload_batch(db, [inp], [])
    assert len(db.select_all('text_ref')) == num_refs,\
        "Duplicate refs allowed to be submitted.."
    return


def test_continuing_upload():
    "Test whether we correcly pick up where we left off when continuing."
    raise SkipTest("Not sure how to actually test this....")


def test_full_remote_upload():
    "Test whether we can perform a targeted upload to aws."
    # This uses a specially curated sample directory designed to access all code
    # paths that the real system might experience, but on a much smaller (thus
    # faster) scale. Errors in the ftp service will not be caught by this test.
    if not TEST_POPULATE:
        raise SkipTest("These feautures are only supported in Python 3.x")
    raise SkipTest("Do not mess with the database...")
    db = remote_startup()
    loc_path = 'test_ftp'
    Medline(ftp_url=loc_path, local=True).populate(db)
    tr_list = db.select_all('text_ref')
    assert all([hasattr(tr, 'pmid') for tr in tr_list]),\
        'All text_refs MUST have pmids by now.'
    PmcOA(ftp_url=loc_path, local=True).populate(db)
    tc_list = db.select_all(
        'text_content',
        db.TextContent.text_type == texttypes.FULLTEXT
        )
    assert len(tc_list), "No fulltext was added."
    Manuscripts(ftp_url=loc_path, local=True).populate(db)
    tc_list = db.select_all(
        'text_content',
        db.TextContent.source == Manuscripts.my_source
        )
    assert len(tc_list), "No manuscripts uploaded."


def test_id_handling_pmc_oa():
    "Test every conceivable combination pmid/pmcid presense."
    if not TEST_POPULATE:
        raise SkipTest("These feautures are only supported in Python 3.x")
    db = local_startup()
    pmc = PmcOA(ftp_url=TEST_FTP, local=True)

    # Initialize with all possible states we could have gotten from medline.
    pm_inp_tpl_list = [
        ('caseA%d' % i, 'PMCcaseA%d' % i) for i in range(2)
        ] + [
        ('caseB%d' % i, None) for i in range(2)
        ] + [
        (None, 'PMCcaseC%d' % i) for i in range(2)
        ] + [
        ('28884161', None),
        ('26977217', 'PMC4771487')
        ]
    db.insert_many(
        'text_ref',
        [dict(zip(('pmid', 'pmcid'), d)) for d in pm_inp_tpl_list]
        )

    # Prepare the 'batch' to be submitted for pmc oa, and try it.
    oa_inp_tpl_list = [
        ('case%s0' % l, 'PMCcase%s0' % l) for l in ['A', 'B', 'C']
        ] + [
        (None, 'PMCcase%s1' % l) for l in ['A', 'B', 'C']
        ] + [
        (None, 'PMC5579538'),  # lookup pmid in db
        (None, 'PMC4238023'),  # lookup no pmid in db
        ('26977217', 'PMC5142709'),  # wrong pmid
        ]
    tr_inp = []
    for pmid, pmcid in oa_inp_tpl_list:
        inp_dict = dict.fromkeys(pmc.tr_cols)
        inp_dict.update(pmcid=pmcid, pmid=pmid)
        tr_inp.append(inp_dict)
    tc_inp = [{'pmcid': pmcid, 'text_type': 'txt', 'content': b'content'}
              for _, pmcid in oa_inp_tpl_list]
    pmc.upload_batch(db, tr_inp, tc_inp)

    # Check the text refs.
    expected_pairs = [
        ('caseA0', 'PMCcaseA0'),
        ('caseA1', 'PMCcaseA1'),
        ('caseB0', 'PMCcaseB0'),
        ('caseB1', None),  # in practice this should be resolved with id_lookup
        ('caseC0', 'PMCcaseC0'),
        (None, 'PMCcaseC1'),
        ('28884161', 'PMC5579538'),
        ('26977217', 'PMC4771487'),
        (None, 'PMCcaseB1'),
        ('25409783', 'PMC4238023')
        ]
    actual_pairs = [(tr.pmid, tr.pmcid) for tr in db.select_all('text_ref')]
    assert_list_equal(actual_pairs, expected_pairs, 'DB text refs incorrect.')

    # Check the text content
    assert len(db.select_all('text_content')) is 8, 'Too much DB text content.'
    return


def test_ftp_service():
    "Test the progenitor childrens' ftp access."
    if not TEST_POPULATE:
        raise SkipTest("These feautures are only supported in Python 3.x")
    cases = [
        ('.csv', 'csv_as_dict'),
        ('.txt', 'file')
        ]
    for Child in [Medline, PmcOA, Manuscripts]:
        c = Child()
        files = c.ftp.ftp_ls()
        for end, meth in cases:
            rel_files = [f for f in files if f.endswith(end)]
            # TODO: check metadata to choose small files.
            if len(rel_files) > 0:
                file = rel_files[0]
                ret = getattr(c.ftp, 'get_' + meth)(file)
                assert ret is not None,\
                    "Failed to load %s from %s." % (file, Child.__name__)
        # TODO: test the download.



