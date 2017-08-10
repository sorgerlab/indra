from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import shutil
import gzip
import lxml.etree.ElementTree as ET
from os import path, walk
from subprocess import call
from collections import namedtuple

from indra.db import insert_text_ref, insert_text_content
from indra.util import zip_string
from indra.tools.reading.pmc_upload_to_s3 import content_type

RE_PATT_TYPE = type(re.compile(''))
# TODO: finish this
SP_INFO = namedtuple('springer_info', ('File', 'date'))

def _find(top_dir, patt):
    '''Find files that match `patt` recursively down from `top_dir`
    
    Note: patt may be a regex string or a regex pattern object.
    '''
    if not isinstance(patt, RE_PATT_TYPE):
        patt = re.compile(patt)

    def match_func(fname):
        return patt.match(fname) is not None

    matches = []
    for root, _, filenames in walk(top_dir):
        for filename in filter(match_func, filenames):
            matches.append(path.join(root, filename))
    return matches


def _pdftotext(pdf_file_path, txt_file_path = None):
    '''Wrapper around the command line function of the same name'''
    if txt_file_path is None:
        txt_file_path = pdf_file_path.replace('.pdf', '.txt')
    elif callable(txt_file_path):
        txt_file_path = txt_file_path(pdf_file_path)

    call(['pdftotext', pdf_file_path, txt_file_path])
    assert path.exists(txt_file_path),\
         "A txt file was not created or name is unknown!"

    return txt_file_path


def upload_springer(springer_dir):
    '''Convert the pdfs to text and upload data to AWS'''
    # TODO: This probably shouldn't be hard-coded.
    txt_path = 'tmp.txt'
    for pdf_path in _find(springer_dir, '.*?\.pdf'):
        pdf_name = path.basename(pdf_path)
        
        # This is where the xml /should/ be.
        art_dirname = path.abspath(pdf_path + '/'.join(4*['..']))
        xml_path_list = _find(art_dirname, pdf_name + '\.xml.*?')
        assert len(xml_path_list) > 0, "We have no metadata"
        if len(xml_path_list) == 1:
            xml = ET.parse(xml_path_list[0])
        elif len(xml_path_list) > 1:
            #TODO: we really should be more intelligent about this
            xml = ET.parse(xml_path_list[0])
        text_ref_id = insert_text_ref(
            source = 'springer',
            doi = xml.find('.//ArticleDOI').text)
        
        #TODO: define the content_type
        content_type = None
        
        txt_path = _pdftotext(pdf_path, txt_path)
        with open(txt_path, 'rb') as f_in:
            with gzip.open(txt_path + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        with open(txt_path, 'rb') as f:
            content = zip_string(f.read())
            
        insert_text_content(text_ref_id, content_type, content)
    return


if __name__ == "__main__":
    #TODO: we should probably support reading from a different
    # directory.
    default_dir = '/groups/lsp/darpa/springer/content/data'
    upload_springer(default_dir)