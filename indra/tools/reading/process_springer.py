from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import shutil
import gzip
from indra.literature import id_lookup, pubmed_client
try:
    import lxml.etree.ElementTree as ET
except:
    import lxml.etree as ET
from os import path, walk, remove
from subprocess import call
from collections import namedtuple

from indra.db import insert_text_ref, insert_text_content
from indra.util import zip_string

RE_PATT_TYPE = type(re.compile(''))
# TODO: finish this
SP_INFO = namedtuple('springer_info', ('File', 'date'))

# TODO: This might do better in util, or somewhere more gnereal ====
def deep_find(top_dir, patt):
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


def pdftotext(pdf_file_path, txt_file_path = None):
    '''Wrapper around the command line function of the same name'''
    if txt_file_path is None:
        txt_file_path = pdf_file_path.replace('.pdf', '.txt')
    elif callable(txt_file_path):
        txt_file_path = txt_file_path(pdf_file_path)

    call(['pdftotext', pdf_file_path, txt_file_path])
    assert path.exists(txt_file_path),\
         "A txt file was not created or name is unknown!"

    return txt_file_path
#====================================================================

def get_xml_data(pdf_path):
    'Get the data from the xml file if present'
    pdf_name = path.basename(pdf_path)
    art_dirname = path.abspath(pdf_path + '/'.join(4*['..']))
    xml_path_list = deep_find(art_dirname, pdf_name.replace('.pdf','\.xml.*?'))
    assert len(xml_path_list) > 0, "We have no metadata"
    if len(xml_path_list) == 1:
        xml = ET.parse(xml_path_list[0])
    elif len(xml_path_list) > 1:
        #TODO: we really should be more intelligent about this
        xml = ET.parse(xml_path_list[0])
    
    # Maybe include the journal subtitle too, in future.
    entry_dict = {'doi':'ArticleDOI',
                  'journal':'JournalTitle',
                  'pub_date':'Year',
                  'publisher':'PublisherName'}
    ref_data = {}
    for table_key, xml_label in entry_dict.items():
        ref_data[table_key] = xml.find('.//' + xml_label).text
    
    return ref_data

def find_other_ids(doi):
    '''Use the doi to try and find the pmid and/or pmcid.'''
    other_ids = dict(zip(['pmid', 'pmcid'], 2*[None]))
    id_dict = id_lookup(doi, 'doi')
    if id_dict['pmid'] is None:
        result_list = pubmed_client.get_ids(doi)
        for res in result_list:
            if 'PMC' in res:
                other_ids['pmcid'] = res
                # This is not very efficient...
                other_ids['pmid'] = id_lookup(res, 'pmcid')['pmid']
                break
        else:
            # This is based on a circumstantial assumption.
            # It will work for the test set, but may fail when
            # upon generalization.
            if len(result_list) == 1:
                other_ids['pmid'] = result_list[0]
            else:
                other_ids['pmid'] = None
    else:
        other_ids['pmid'] = id_dict['pmid']
        other_ids['pmcid'] = id_dict['pmcid']
        
    return other_ids

def process_one_pdf(pdf_path, txt_path):
    'Convert the pdf to txt and zip it'
    txt_path = pdftotext(pdf_path, txt_path)
    with open(txt_path, 'rb') as f_in:
        with gzip.open(txt_path + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    with open(txt_path, 'rb') as f:
        content = zip_string(f.read().decode('utf-8'))
    remove(txt_path) # Only a tmp file.
    return content


def upload_springer(springer_dir):
    '''Convert the pdfs to text and upload data to AWS
    
    Note: Currently does nothing.
    '''
    # TODO: We should probably filter which articles we process
    txt_path = 'tmp.txt'
    uploaded = []
    for pdf_path in deep_find(springer_dir, '.*?\.pdf'):
        ref_data = get_xml_data(pdf_path)
        ref_data.update(find_other_ids(ref_data['doi']))
        
        if ref_data['pmid'] is None and ref_data['pmcid'] is None:
            # We will for now assume this article is not relevant.
            continue
                
        #text_ref_id = insert_text_ref(source = 'springer', **ref_data)
        content_type = None #TODO: define the content_type
        content = process_one_pdf(pdf_path, txt_path)
        #insert_text_content(text_ref_id, content_type, content)
        uploaded.append(pdf_path)
    return uploaded


if __name__ == "__main__":
    #TODO: we should probably support reading from a different
    # directory.
    default_dir = '/groups/lsp/darpa/springer/content/data'
    upload_springer(default_dir)