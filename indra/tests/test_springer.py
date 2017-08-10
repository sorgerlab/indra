from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.tools.reading import process_springer as ps

TOP_DIR = 'springer_mock'
PDF_PATH = TOP_DIR + '/test_article_1v/BodyRef/PDF/15010_2002_Article_1083.pdf'

def test_xml_read():
    'Tests whether the xml files are being read'
    ps.get_xml_data(PDF_PATH)
    return
    
def test_convert_and_zip():
    'Tests that the pdf can be converted and zipped'
    ps.process_one_pdf(PDF_PATH, 'tmp_test.txt')
    return
    
def test_upload():
    'Tests whether the file was uploaded successfully'
    ps.upload_springer(TOP_DIR)
    return