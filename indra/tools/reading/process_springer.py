from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import shutil
import gzip
import lxml.etree.ElementTree as ET
from os import path, walk, remove
from subprocess import call
from collections import namedtuple

from indra.db import insert_text_ref, insert_text_content
from indra.util import zip_string

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

def get_xml_data(pdf_path):
    'Get the data from the xml file if present'
    pdf_name = path.basename(pdf_path)
    art_dirname = path.abspath(pdf_path + '/'.join(4*['..']))
    xml_path_list = _find(art_dirname, pdf_name + '\.xml.*?')
    assert len(xml_path_list) > 0, "We have no metadata"
    if len(xml_path_list) == 1:
        xml = ET.parse(xml_path_list[0])
    elif len(xml_path_list) > 1:
        #TODO: we really should be more intelligent about this
        xml = ET.parse(xml_path_list[0])
    
    entry_dict = {'doi':'ArticleDOI',
                  'journal':['JournalTitle', 'JournalSubTitle'],
                  'pub_date':'Year',
                  'publisher':'PublisherName'}
    ref_data = {}
    for table_key, xml_label in entry_dict:
        ref_data[table_key] = xml.find('.//' + xml_label)
    
    return ref_data

def process_one_pdf(pdf_path, txt_path):
    'Convert the pdf to txt and zip it'
    txt_path = _pdftotext(pdf_path, txt_path)
    with open(txt_path, 'rb') as f_in:
        with gzip.open(txt_path + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    with open(txt_path, 'rb') as f:
        content = zip_string(f.read())
    remove(txt_path) # Only a tmp file.
    return content


def upload_springer(springer_dir):
    '''Convert the pdfs to text and upload data to AWS'''
    # TODO: We should probably filter which articles we process
    txt_path = 'tmp.txt'
    for pdf_path in _find(springer_dir, '.*?\.pdf'):
        ref_data = get_xml_data(pdf_path)
        text_ref_id = insert_text_ref(source = 'springer', **ref_data)
        content_type = None #TODO: define the content_type
        content = process_one_pdf(pdf_path, txt_path)
        insert_text_content(text_ref_id, content_type, content)
    return


if __name__ == "__main__":
    #TODO: we should probably support reading from a different
    # directory.
    default_dir = '/groups/lsp/darpa/springer/content/data'
    upload_springer(default_dir)