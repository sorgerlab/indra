from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import re

from indra.tools.reading import process_springer as ps
from os.path import exists, dirname, join, realpath, basename
from time import sleep
from datetime import datetime
from os import makedirs, removedirs, remove
from nose import SkipTest
from indra.db import DatabaseManager

TOP_DIR = join(realpath(dirname(__file__)), 'springer_mock')
if not exists(TOP_DIR):
    raise SkipTest("Cannot run these tests without mock directory.")

PDF_PATH = TOP_DIR + '/ART_1_NOPMID/BodyRef/PDF/15010_2002_Article_1083.pdf'
FINDME = '{}.findme.tmp'

TEST_DB_FILE = 'springer_test.db'
TEST_DB_HOST = 'sqlite:///' + TEST_DB_FILE


def test_basic_deep_find():
    'Test the basic functionality of the deep_find method'
    # Setup
    findables = [
        FINDME.format('shallow'),
        join('inhere', FINDME.format('deep'))
        ]
    try:
        for findable in findables:
            find_dir = dirname(findable)
            if find_dir != '' and not exists(find_dir):
                makedirs(find_dir)
            open(findable, 'w').close()

        # Test
        flist = ps.deep_find('.', FINDME.format('.*?'))
        assert set(map(basename, flist)) == set(map(basename, findables)),\
            "Expect to find: %s, but found %s." % (str(findables), str(flist))
    finally:
        for findable in findables:
            if exists(findable):
                remove(findable)
            find_dif = dirname(findable)
            if find_dif != '' and exists(find_dif):
                removedirs(find_dif)
    return


def test_time_sensitive_deep_find():
    'Tests whether the only recent files can be searched.'
    old_fname = join(FINDME.format('old'))
    open(old_fname, 'w').close()
    sleep(1)
    check_from_date = datetime.now()
    sleep(1)
    new_fname = FINDME.format('new')
    open(new_fname, 'w').close()
    try:
        flist = ps.deep_find(
            '.', 
            FINDME.format('.*?'),
            since_date=check_from_date
            )
        assert len(flist) == 1 and new_fname in set(map(basename, flist)),\
            'Not just recent files found: %s.' % str(flist)
    finally:
        for fname in [old_fname, new_fname]:
            if exists(fname):
                remove(fname)
    return


def test_xml_read():
    'Tests whether the xml files are being read'
    ps.get_xml_data(PDF_PATH, entry_dict={'ref_data': {'doi': 'ArticleDOI'}})
    return


def test_convert_and_zip():
    'Tests that the pdf can be converted and zipped'
    ps.process_one_pdf(PDF_PATH, 'tmp_test.txt')
    assert not exists('tmp_test.txt'), 'Temp file not so temporary!'
    return


def test_id_finding():
    'Tests the `find_other_ids` function'
    def test_a_doi(doi, expect_ids, err_str):
        test_ids = ps.find_other_ids(doi)
        assert test_ids == expect_ids, err_str

    # TODO: Test a doi that returns no other ids
    # this requires that I find one which satisfies.
    test_a_doi(
        '10.2478/s11686-011-0061-7',
        dict(pmid=None, pmcid=None),
        'Should NOT find a pmid OR a pmcid.'
        )

    # Test a doi that has only a pmid
    # FIXME: Find an actual example. The pmid found was actually erroneous.
    # test_a_doi(
    #    '10.1007/s15010-002-1083-8',
    #    dict(pmid = '26811110', pmcid = None),
    #    'Should ONLY find a pmid.'
    #    )

    # Test a doi that has a pmid and a pmcid
    test_a_doi(
        '10.1007/s11606-009-1182-7',
        dict(pmid='19967464', pmcid='PMC2839328'),
        'Should find a pmid AND a pmcid.'
        )
    return


def test_upload():
    'Tests whether the file was uploaded successfully'
    # Prepare the tables.
    # TODO: Do some setup in the database, to test whether we handle
    # pre-existing content correctly.
    db = DatabaseManager(TEST_DB_HOST, sqltype='sqlite')
    if not len(db.get_active_tables()):
        db.create_tables()
        created_tables = True
    else:
        db._clear()
        created_tables = False

    try:
        # Upload the test documents.
        uploads = ps.upload_springer(
            TOP_DIR,
            host=TEST_DB_HOST,
            sqltype='sqlite',
            verbose=True)
        art_patt = re.compile('.*?(ART.*?).BodyRef.*?')
        uploads_arts = set()
        for upload in uploads:
            uploads_arts.add(art_patt.match(upload['path']).groups()[0])
        expected_uploads = {'ART_PMCID', 'ART_PMID'}
        assert uploads_arts == expected_uploads,\
            'Wrong articles uploaded.\n' +\
            'Expected: %s,\ngot: %s.' % (str(expected_uploads), str(uploads_arts))

        # Check that the correct database entries were made.
        db.grab_session()
        tref_dict = {}
        tcont_dict = {}
        for upload in uploads:
            key = basename(upload['path'])
            tref = db.select_one('text_ref', db.TextRef.doi == upload['doi'])
            tref_dict[key] = tref
            if tref is None:
                continue  # Cannot look for content w/o tref.id
            tcont = db.select('text_content', db.TextContent.text_ref == tref)
            tcont_dict[key] = tcont if len(tcont) > 0 else None
        err_fmt = "Did not insert all text %s.\nGot some \'None\'s: %s."
        filter_dict = lambda d: filter(lambda _,v: v is None, d.items())
        assert None not in list(tref_dict.values()),\
            err_fmt % ('ref', filter_dict(tref_dict))
        assert None not in list(tcont_dict.values()),\
            err_fmt % ('content', filter_dict(tcont_dict)) 

    finally:
        if created_tables:
            db.drop_tables()
        else:
            db._clear()
    return
