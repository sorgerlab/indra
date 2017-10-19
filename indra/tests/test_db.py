from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import re
from os import remove, path
from sys import version_info
from nose import SkipTest
from nose.tools import assert_equal
from functools import wraps
from sqlalchemy.exc import IntegrityError
from indra.db import DatabaseManager, texttypes, get_defaults
IS_PY3 = True
if version_info.major is not 3:
    IS_PY3 = False
if IS_PY3:
    from indra.db.manage_content import Medline, PmcOA, Manuscripts

defaults = get_defaults()
test_defaults = {k: v for k, v in defaults.items() if 'test' in k}

# Get host for the test database from system defaults.
# TODO: implement setup-teardown system.
TEST_HOST = None
TEST_HOST_TYPE = ''
key_list = list(test_defaults.keys())
key_list.sort()
for k in key_list:
    v = test_defaults[k]
    m = re.match('(\w+)://.*?/([\w.]+)', v)
    if m is None:
        raise SkipTest(
            "No hosts found for testing. Please update defaults.txt if you\n"
            "wish to run these tests. Test host entries must be of the form\n"
            "test*=<sql type>://<passwords and such>/<database name>")
    sqltype = m.groups()[0]
    db_name = m.groups()[1]
    try:
        db = DatabaseManager(v, sqltype=sqltype)
        db._clear()
    except Exception as e:
        print("Tried to use %s, but failed due to:\n%s" % (k, e))
        continue  # Clearly this test table won't work.
    if db_name.endswith('.db'):
        TEST_FILE = db_name
        if path.exists(TEST_FILE):
            remove(TEST_FILE)
    else:
        TEST_FILE = None
    TEST_HOST = v
    TEST_HOST_TYPE = sqltype
    print("Using test database %s." % k)
    break
else:
    raise SkipTest("Not able to start up any of the available test hosts.")

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


def get_db(clear=True):
    "Set up the database for testing."
    db = DatabaseManager(TEST_HOST, sqltype=TEST_HOST_TYPE)
    db.grab_session()
    if clear:
        db._clear()
    return db


def needs_py3(func):
    @wraps(func)
    def test_with_py3_func(*args, **kwargs):
        if not IS_PY3:
            raise SkipTest("This tests feautures only supported in Python 3.x")
        return func(*args, **kwargs)
    return test_with_py3_func


#==============================================================================
# The following are tests for the database manager itself.
#==============================================================================
def test_create_tables():
    "Test the create_tables feature"
    db = get_db()
    db.create_tables()
    assert_contents_equal(db.get_active_tables(), db.get_tables())


def test_insert_and_query_pmid():
    "Test that we can add a text_ref and get the text_ref back."
    db = get_db()
    pmid = '1234'
    text_ref_id = db.insert('text_ref', pmid=pmid)
    entries = db.select_all('text_ref', db.TextRef.pmid == pmid)
    assert_equal(len(entries), 1, "One item inserted, multiple entries found.")
    assert_equal(entries[0].pmid, pmid)
    assert_equal(entries[0].id, text_ref_id, "Got back wrong text_ref_id.")


def test_uniqueness_text_ref_doi_pmid():
    "Test uniqueness enforcement behavior for text_ref insertion."
    db = get_db()
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
    db = get_db()
    url = 'http://foobar.com'
    db.insert('text_ref', url=url)
    try:
        db.insert('text_ref', url=url)
    except IntegrityError:
        return  # PASS
    assert False, "Uniqueness was not enforced."


def test_get_abstracts():
    "Test the ability to get a list of abstracts."
    db = get_db()

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
    db = get_db()
    db.insert_many('text_ref', [{'pmid': '1234'}, {'pmid': '5678'}])
    pmid_list = db.get_all_pmids()
    assert_contents_equal(pmid_list, ['1234', '5678'])


#==============================================================================
# The following are tests for the initial population of the database. This
# includes uploading data from Medline, PMC, Springer, and Elsevier. These
# tend to make greater use of the database, and are likely to be slower.
#==============================================================================
TEST_FTP = path.abspath(path.join(path.dirname(__file__), 'test_ftp'))
if IS_PY3 and not path.exists(TEST_FTP):
    print("Creating test directory. This could take a while...")
    from indra.db.build_sample_set import build_set
    build_set(2, TEST_FTP)


@needs_py3
def test_full_upload():
    "Test whether we can perform a targeted upload to a test db."
    # This uses a specially curated sample directory designed to access most
    # code paths that the real system might experience, but on a much smaller
    # (thus faster) scale. Errors in the ftp service will not be caught by
    # this test.
    db = get_db()
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
    tc_list = db.select_all('text_content')
    set_exp = {('manuscripts', 'xml', 'fulltext'),
               ('pmc_oa', 'xml', 'fulltext'),
               ('pubmed', 'text', 'abstract')}
    set_got = set([(tc.source, tc.format, tc.text_type) for tc in tc_list])
    assert set_exp == set_got,\
        "Expected %s, got %s for content layout." % (set_exp, set_got)


@needs_py3
def test_multiple_pmids():
    "Test that pre-existing pmids are correctly handled."
    db = get_db()
    med = Medline(ftp_url=TEST_FTP, local=True)
    med.populate(db)
    num_refs = len(db.select_all('text_ref'))
    med.populate(db)
    assert len(db.select_all('text_ref')) == num_refs,\
        "Duplicate pmids allowed to be submitted.."
    return


@needs_py3
def test_multible_pmc_oa_content():
    "Test to make sure repeated content is handled correctly."
    db = get_db()
    pmc = PmcOA(ftp_url=TEST_FTP, local=True)
    pmc.populate(db)
    num_conts = len(db.select_all('text_content'))
    pmc.populate(db)
    assert len(db.select_all('text_content')) == num_conts,\
        "Duplicate text content allowed to be submitted."
    return


@needs_py3
def test_multiple_text_ref_pmc_oa():
    "Test whether a duplicate text ref in pmc oa is handled correctly."
    db = get_db()
    pmc = PmcOA(ftp_url=TEST_FTP, local=True)
    inp = dict.fromkeys(pmc.tr_cols)
    inp.update(pmcid='PMC5579538', doi='10.1021/acsomega.7b00205')
    pmc.upload_batch(db, [inp], [])
    num_refs = len(db.select_all('text_ref'))
    pmc.upload_batch(db, [inp], [])
    assert len(db.select_all('text_ref')) == num_refs,\
        "Duplicate refs allowed to be submitted.."
    return


@needs_py3
def test_id_handling_pmc_oa():
    "Test every conceivable combination pmid/pmcid presense."
    db = get_db()
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
    assert_contents_equal(
        actual_pairs,
        expected_pairs,
        'DB text refs incorrect.'
        )

    # Check the text content
    assert len(db.select_all('text_content')) is 8, 'Too much DB text content.'
    return


@needs_py3
def test_ftp_service():
    "Test the NIH FTP access client on the content managers."
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
                fname = rel_files[0]
                ret = getattr(c.ftp, 'get_' + meth)(fname)
                assert ret is not None,\
                    "Failed to load %s from %s." % (fname, Child.__name__)
        # TODO: test the download.
