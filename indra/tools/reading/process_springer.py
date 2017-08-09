from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import re

from os import listdir, path, walk
from subprocess import call

RE_PATT_TYPE = type(re.compile(''))

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


def find(top_dir, patt):
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


def convert_pdfs(dir_path, name_rule):
    '''Convert the pdfs in the springer directory to text.'''
    pdf_path_list = find(dir_path, '.*?\.pdf')
    txt_path_list = []
    for pdf_path in pdf_path_list:
        print("Converting %s to txt." % path.basename(pdf_path))
        txt_path_list = pdftotext(pdf_path, name_rule)
    return txt_path_list

def default_name_rule(pdf_path):
    '''This is a default name rule for use with convert_pdfs'''
    pass # TODO: Actually write this. Must have access to group data first.


def add_to_db():
    '''Add the text files to the database'''
    # Use XML files for metadata
    pass