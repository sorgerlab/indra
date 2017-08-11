from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import re

from indra.tools.reading import process_springer as ps
from os.path import exists

TOP_DIR = 'springer_mock'
PDF_PATH = TOP_DIR + '/ART_1_NOPMID/BodyRef/PDF/15010_2002_Article_1083.pdf'

def test_xml_read():
    'Tests whether the xml files are being read'
    ps.get_xml_data(PDF_PATH)
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