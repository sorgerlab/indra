from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import re
from os import remove, path
from sys import version_info, argv
from nose import SkipTest
from nose.tools import assert_equal
from functools import wraps
from sqlalchemy.exc import IntegrityError
from indra.db.database_manager import DatabaseManager
from indra.db.util import get_abstracts_by_pmids, get_defaults, NestedDict
from nose.plugins.attrib import attr
from indra.db.reading_manager import BulkReadingManager, BulkLocalReadingManager

IS_PY3 = True
if version_info.major is not 3:
    IS_PY3 = False
if IS_PY3:
    from indra.db.content_manager import Pubmed, PmcOA, Manuscripts, Elsevier

if '-a' in argv:
    attr_str = argv[argv.index('-a')+1]
    if any([not_attr in attr_str for not_attr in
            ('!nonpublic', '!webservice')]):
        raise SkipTest("Every test is nonpublic and a webservice.")

defaults = get_defaults()
test_defaults = {k: v for k, v in defaults.items() if 'test' in k}

# Get host for the test database from system defaults.
TEST_HOST = None
TEST_HOST_TYPE = ''
key_list = list(test_defaults.keys())
key_list.sort()
report = ''
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
        db._clear(force=True)
    except Exception as e:
        report += "Tried to use %s, but failed due to:\n%s\n" % (k, e)
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
    raise SkipTest("Not able to start up any of the available test hosts:\n"
                   + report)


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


def capitalize_list_of_tpls(l):
    return [tuple([i.upper() if isinstance(i, str) else i for i in e])
            for e in l]


def get_db(clear=True):
    "Set up the database for testing."
    db = DatabaseManager(TEST_HOST, sqltype=TEST_HOST_TYPE)
    db.grab_session()
    if clear:
        db._clear(force=True)
    return db


def needs_py3(func):
    @wraps(func)
    def test_with_py3_func(*args, **kwargs):
        if not IS_PY3:
            raise SkipTest("This tests features only supported in Python 3.x")
        return func(*args, **kwargs)
    return test_with_py3_func


@needs_py3
def get_db_with_pubmed_content():
    "Populate the database with sample content from pubmed."
    db = get_db()
    Pubmed(ftp_url=TEST_FTP, local=True).populate(db)
    return db


@needs_py3
def get_db_with_ftp_content():
    "Populate database with content from all the ftp services"
    db = get_db_with_pubmed_content()
    PmcOA(ftp_url=TEST_FTP, local=True).populate(db)
    Manuscripts(ftp_url=TEST_FTP, local=True).populate(db)
    return db


#==============================================================================
# The following are tests for the database manager itself.
#==============================================================================
@attr('nonpublic')
def test_create_tables():
    "Test the create_tables feature"
    db = get_db()
    db.create_tables()
    assert_contents_equal(db.get_active_tables(), db.get_tables())


@attr('nonpublic')
def test_insert_and_query_pmid():
    "Test that we can add a text_ref and get the text_ref back."
    db = get_db()
    pmid = '1234'
    text_ref_id = db.insert('text_ref', pmid=pmid)
    entries = db.select_all('text_ref', db.TextRef.pmid == pmid)
    assert_equal(len(entries), 1, "One item inserted, multiple entries found.")
    assert_equal(entries[0].pmid, pmid)
    assert_equal(entries[0].id, text_ref_id, "Got back wrong text_ref_id.")


@attr('nonpublic')
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
        db._clear(force=True)
    assert False, "Uniqueness was not enforced."


@attr('nonpublic')
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


@attr('nonpublic')
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
    received = get_abstracts_by_pmids(db, ['1234', '5678', '1357'],
                                      unzip=False)
    assert_contents_equal(expected, received, "Didn't get expected abstracts.")


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
@attr('nonpublic', 'slow')
def test_full_upload():
    "Test whether we can perform a targeted upload to a test db."
    # This uses a specially curated sample directory designed to access most
    # code paths that the real system might experience, but on a much smaller
    # (thus faster) scale. Errors in the ftp service will not be caught by
    # this test.

    # Test the medline/pubmed upload.
    db = get_db_with_pubmed_content()
    tr_list = db.select_all('text_ref')
    assert len(tr_list), "No text refs were added..."
    assert all([hasattr(tr, 'pmid') for tr in tr_list]),\
        'All text_refs MUST have pmids by now.'

    # Test the pmc oa upload.
    PmcOA(ftp_url=TEST_FTP, local=True).populate(db)
    tcs_pmc = db.filter_query(
        db.TextContent,
        db.TextContent.source == PmcOA.my_source).count()
    assert tcs_pmc, "No pmc oa fulltext was added."
    trs_w_pmcids = db.filter_query(
        db.TextRef,
        db.TextRef.pmcid.isnot(None)).count()
    assert trs_w_pmcids >= tcs_pmc,\
        "Only %d of at least %d pmcids added." % (trs_w_pmcids, tcs_pmc)

    # Test the manuscripts upload.
    Manuscripts(ftp_url=TEST_FTP, local=True).populate(db)
    tcs_manu = db.filter_query(
        db.TextContent,
        db.TextContent.source == Manuscripts.my_source
        ).count()
    assert tcs_manu, "No manuscripts uploaded."
    trs_w_mids = db.filter_query(
        db.TextRef,
        db.TextRef.manuscript_id.isnot(None)
        ).count()
    assert trs_w_mids >= tcs_manu,\
        "Only %d of at least %d manuscript ids added." % (trs_w_mids, tcs_manu)

    # Some overal checks.
    tc_list = db.select_all(db.TextContent)
    set_exp = {('manuscripts', 'xml', 'fulltext'),
               ('pmc_oa', 'xml', 'fulltext'),
               ('pubmed', 'text', 'abstract')}
    set_got = set([(tc.source, tc.format, tc.text_type) for tc in tc_list])
    assert set_exp == set_got,\
        "Expected %s, got %s for content layout." % (set_exp, set_got)

    # Test careful upload of medline (very shallow test...checks only for
    # critical failures)
    m = Pubmed(ftp_url=TEST_FTP, local=True)
    m.load_files(db, 'baseline', carefully=True)


@needs_py3
@attr('nonpublic')
def test_multiple_pmids():
    "Test that pre-existing pmids are correctly handled."
    db = get_db()
    med = Pubmed(ftp_url=TEST_FTP, local=True)
    med.populate(db)
    num_refs = len(db.select_all('text_ref'))
    med.populate(db)
    assert len(db.select_all('text_ref')) == num_refs,\
        "Duplicate pmids allowed to be submitted.."
    return


@needs_py3
@attr('nonpublic')
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
@attr('nonpublic')
def test_multiple_text_ref_pmc_oa():
    "Test whether a duplicate text ref in pmc oa is handled correctly."
    db = get_db()
    pmc = PmcOA(ftp_url=TEST_FTP, local=True)
    pmc.review_fname = 'test_review_multiple_text_ref_pmc_oa.txt'
    inp = dict.fromkeys(pmc.tr_cols)
    inp.update(pmcid='PMC5579538', doi='10.1021/acsomega.7b00205')
    pmc.upload_batch(db, [inp], [])
    num_refs = len(db.select_all('text_ref'))
    pmc.upload_batch(db, [inp], [])
    assert len(db.select_all('text_ref')) == num_refs,\
        "Duplicate refs allowed to be submitted.."
    remove(pmc.review_fname)
    return


@needs_py3
@attr('nonpublic')
def test_id_handling_pmc_oa():
    "Test every conceivable combination pmid/pmcid presence."
    db = get_db()
    pmc = PmcOA(ftp_url=TEST_FTP, local=True)

    # Initialize with all possible states we could have gotten from medline.
    pm_inp_tpl_list = capitalize_list_of_tpls([
        ('caseA%d' % i, 'PMCcaseA%d' % i) for i in range(2)
        ] + [
        ('caseB%d' % i, None) for i in range(2)
        ] + [
        (None, 'PMCcaseC%d' % i) for i in range(2)
        ] + [
        ('caseMisMatchA', 'PMCcaseMisMatchB'),
        ('caseMisMatchB', 'PMCcaseMisiMatchB'),
        ('caseMultiMatch', 'PMCcaseMultiMatch'),
        ('28884161', None),
        ('26977217', 'PMC4771487')
        ])
    db.insert_many(
        'text_ref',
        [dict(zip(('pmid', 'pmcid'), d)) for d in pm_inp_tpl_list]
        )

    # Prepare the 'batch' to be submitted for pmc oa, and try it.
    oa_inp_tpl_list = capitalize_list_of_tpls([
        ('case%s0' % l, 'PMCcase%s0' % l) for l in ['A', 'B', 'C']
        ] + [
        (None, 'PMCcase%s1' % l) for l in ['A', 'B', 'C']
        ] + [
        (None, 'PMC5579538'),  # lookup pmid in db
        (None, 'PMC4238023'),  # lookup no pmid in db
        ('26977217', 'PMC5142709'),  # conflicting pmcid
        ('caseMisMatchB', 'PMCcaseMisMatchA'),  # multiple matches
        ('caseMultiMatch', 'PMCnotmatching'),
        ('notmatching', 'PMCcaseMultiMatch'),
        ])
    tr_inp = []
    for pmid, pmcid in oa_inp_tpl_list:
        inp_dict = dict.fromkeys(pmc.tr_cols)
        inp_dict.update(pmcid=pmcid, pmid=pmid)
        tr_inp.append(inp_dict)
    tc_inp = [{'pmcid': pmcid, 'text_type': 'txt', 'content': b'content'}
              for _, pmcid in oa_inp_tpl_list]
    pmc.review_fname = 'test_review_%s.txt' % pmc.my_source
    pmc.upload_batch(db, tr_inp, tc_inp)

    # Check the text refs.
    expected_pairs = capitalize_list_of_tpls([
        ('caseA0', 'PMCcaseA0'),
        ('caseA1', 'PMCcaseA1'),
        ('caseB0', 'PMCcaseB0'),
        ('caseB1', None),  # in practice this should be resolved with id_lookup
        ('caseC0', 'PMCcaseC0'),
        (None, 'PMCcaseC1'),
        ('28884161', 'PMC5579538'),
        ('26977217', 'PMC4771487'),
        (None, 'PMCcaseB1'),
        ('25409783', 'PMC4238023'),
        ('caseMisMatchA', 'PMCcaseMisMatchB'),
        ('caseMisMatchB', 'PMCcaseMisiMatchB'),
        ('caseMultiMatch', 'PMCcaseMultiMatch'),
        ])
    actual_pairs = [(tr.pmid, tr.pmcid) for tr in db.select_all('text_ref')]
    assert_contents_equal(actual_pairs, expected_pairs,
                          'DB text refs incorrect.')

    with open(pmc.review_fname, 'r') as f:
        found_conflict_msg = False
        for line in f.read().splitlines():
            if all([word in line for word in
                    ['PMC4771487', 'PMC5142709', 'conflicting pmcid']]):
                found_conflict_msg = True
                break
        assert found_conflict_msg

    # Check the text content
    assert len(db.select_all('text_content')) is 8, 'Too much DB text content.'
    remove(pmc.review_fname)
    return


@needs_py3
@attr('nonpublic')
def test_medline_ref_checks():
    "Test the text ref checks used by medline."
    db = get_db()
    med = Pubmed(ftp_url=TEST_FTP, local=True)

    def check_input(input_pairs, expected_pairs, carefully, num):
        article_info = {pmid: dict(zip(['pmid', 'pmcid'], [pmid, pmcid]))
                        for pmid, pmcid in input_pairs}
        med.load_text_refs(db, article_info, carefully)
        actual_pairs = [(tr.pmid, tr.pmcid)
                        for tr in db.select_all(db.TextRef)]
        desc = 'careful' if carefully else 'careless'
        msg = 'DB text refs mismatch after upload %d (%s)' % (num, desc)
        actual_pairs.sort(key=str)
        expected_pairs.sort(key=str)
        assert_contents_equal(actual_pairs, expected_pairs, msg)

    expected_pairs = [
        ('CASEA', None),
        ('CASEB', 'PMCIDCASEB'),
        ('CASEC', None),
        ('CASED', 'PMCIDCASED')
        ]

    # Upload round 1
    check_input(
        [
            ('CASEA', None),
            ('CASEB', 'PMCIDCASEB'),
            ('CASEC', None),
            ('CASEC', None),
            ('CASED', None),
            ('CASED', 'PMCIDCASED')
            ],
        expected_pairs,
        False,
        1
        )

    # Upload round 2
    expected_pairs += [
        ('CASEE', None)
        ]
    check_input(
        [
            ('CASEE', None),
            ('CASEC', 'PMCIDCASEC'),
            ('CASEH1', 'PMCIDCASEH'),
            ('CASEK', 'PMCIDCASEK1')
            ],
        expected_pairs + [
            ('CASEH1', 'PMCIDCASEH'),
            ('CASEK', 'PMCIDCASEK1')
            ],
        False,
        2
        )

    # Interlude
    db.insert_many('text_ref', [
        {'pmcid': 'PMCIDCASEG'},
        ])

    # Upload round 3
    input_pairs = expected_pairs + [
        ('CASEF', None),
        ('CASEC', 'PMCIDCASEC'),
        ('CASEG', 'PMCIDCASEG'),
        ('CASEH2', 'PMCIDCASEH'),  # this should trigger a review.
        ('CASEK', 'PMCIDCASEK2')  # and so should this
        ]
    expected_pairs.remove(('CASEC', None))
    expected_pairs += [
        ('CASEF', None),
        ('CASEC', 'PMCIDCASEC'),
        ('CASEG', 'PMCIDCASEG'),
        ('CASEH1', 'PMCIDCASEH'),
        ('CASEK', 'PMCIDCASEK1')
        ]
    med.review_fname = 'test_review_%s.txt' % med.my_source
    open(med.review_fname, 'a+').close()
    with open(med.review_fname, 'r') as f:
        num_orig_lines = len(f.readlines())
    check_input(
        input_pairs,
        expected_pairs,
        True,
        3
        )
    with open(med.review_fname, 'r') as f:
        lines = f.readlines()
        assert len(lines) == num_orig_lines + 2, \
            "Not all new reviews added: %d / %d" % (len(lines),
                                                    num_orig_lines + 2)
    remove(med.review_fname)
    return


@needs_py3
@attr('nonpublic')
def test_elsevier_upload():
    "Test that we can upload elsevier content."
    db = get_db_with_ftp_content()
    Elsevier().populate(db)
    up_q = db.filter_query(
        db.Updates,
        db.Updates.source == Elsevier.my_source
        )
    num_updates = up_q.count()
    assert num_updates == 1, "Got %d updates, not 1." % num_updates
    assert up_q.all()[0].init_upload, \
        "Update entry not listed as initial upload."
    tc_q = db.filter_query(
        db.TextContent,
        db.TextContent.source == Elsevier.my_source
        )
    num_elsevier = tc_q.count()
    assert num_elsevier > 0, "Got no elsevier content."


@needs_py3
@attr('nonpublic', 'slow')
def test_sparser_initial_reading():
    "Test the initial reading of of sparser content"
    db = get_db_with_ftp_content()
    BulkLocalReadingManager('sparser', n_proc=1).read_all(db)
    sparser_updates_q = db.filter_query(db.ReadingUpdates,
                                        db.ReadingUpdates.reader == 'SPARSER')
    assert sparser_updates_q.count() == 1, "Update was not logged."
    sparser_readings_q = db.filter_query(db.Readings,
                                         db.Readings.reader == 'SPARSER')
    assert sparser_readings_q.count() > 0, "Failed to produce readings."
    sparser_stmts_q = db.filter_query(db.Statements,
                                      db.Statements.reader_ref == db.Readings.id,
                                      db.Readings.reader == 'SPARSER')
    assert sparser_stmts_q.count() > 0


def test_nested_dict():
    d = NestedDict()
    print(d)
    d['A']['B']['C'] = 3
    d['B']['C'] = {'D': 2}
    print(d)
    assert d['A']['B']['C'] == 3
    assert d.get('A') == d['A']
    assert d.gets('A') == [d['A']]
    assert d.get('C') in [3, {'D': 2}]
    assert d.get('D') == 2
    assert d.get_path('C') in [(('A', 'B', 'C'), 3), (('B', 'C'), {'D': 2})]
    assert_contents_equal([str(v) for v in d.gets('C')],
                          ['3', str(d['B']['C'])])
    d.export_dict()  # Should probably test for matching contents
    assert_contents_equal([str(v) for v in d.get_paths('C')],
                          [str((('A', 'B', 'C'), 3)),
                           str((('B', 'C'), d['B']['C']))])


@needs_py3
@attr('nonpublic', 'slow')
def test_ftp_service():
    "Test the NIH FTP access client on the content managers."
    cases = [
        ('.csv', 'csv_as_dict'),
        ('.txt', 'file')
        ]
    for Child in [Pubmed, PmcOA, Manuscripts]:
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
