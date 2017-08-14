from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import re

from indra.tools.reading import process_springer as ps
from os.path import exists, dirname, join, realpath, basename
from time import sleep
from datetime import datetime
from os import makedirs, removedirs, remove

TOP_DIR = join(realpath(dirname(__file__)), 'springer_mock')
PDF_PATH = TOP_DIR + '/ART_1_NOPMID/BodyRef/PDF/15010_2002_Article_1083.pdf'
FINDME = '{}.findme.tmp'

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
        def clean(fname_list):
            return [basename(f) for f in fname_list]
        assert clean(flist) == clean(findables),\
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
    old_fname = FINDME.format('old')
    open(old_fname, 'w').close()
    sleep(5)
    check_from_date = datetime.now()
    sleep(1)
    new_fname = FINDME.format('old')
    open(new_fname, 'w').close()
    try:
        assert exists(old_fname) and exists(new_fname),\
            'Files to be found not created.'
        flist = ps.deep_find(
            '.', 
            FINDME.format('.*?'), 
            since_date = check_from_date
            )
        assert len(flist) == 1 and new_fname in flist,\
            'Not just recent files found: %s.' % str(flist)
    finally:
        for fname in [old_fname, new_fname]:
            if exists(fname):
                remove(fname)
    return

def test_xml_read():
    'Tests whether the xml files are being read'
    ps.get_xml_data(PDF_PATH, entry_dict = {'ref_data':{'doi':'ArticleDOI'}})
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
    test_a_doi(
        '10.1007/s15010-002-1083-8',
        dict(pmid = '26811110', pmcid = None),
        'Should ONLY find a pmid.'
        )
        
    # Test a doi that has a pmid and a pmcid
    test_a_doi(
        '10.1007/s11606-009-1182-7',
        dict(pmid='19967464', pmcid='PMC2839328'),
        'Should find a pmid AND a pmcid.'
        )
    return
    
def test_upload():
    'Tests whether the file was uploaded successfully'
    uploads = ps.upload_springer(TOP_DIR)
    art_patt = re.compile('.*?(ART.*?).BodyRef.*?')
    uploads_arts = set()
    for upload in uploads:
        uploads_arts.add(art_patt.match(upload).groups()[0])
    expected_uploads = {'ART_1_NOPMID', 'ART_2_NOPMID', 'ART_PMCID', 'ART_PMID'}
    assert uploads_arts ==  expected_uploads,\
        'Wrong articles uploaded.' +\
        'Expected: %s, got: %s.'%(str(expected_uploads), str(uploads_arts)) 
    return