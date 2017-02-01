from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import glob
from indra import sparser

base_folder = os.path.join(os.environ['HOME'],
                           'data/darpa/phase3_eval/sources/sparser-20170127')

def get_file_names(base_dir):
    fnames = glob.glob(os.path.join(base_dir, '*.xml'))
    return fnames

def get_file_stmts(fname):
    with open(fname, 'rt') as fh:
        xml_str = fh.read()
        xml_str = xml_str.replace('hmsid', 'pmid')
        sp = sparser.process_xml(xml_str)
        if sp is None:
            return []
        return sp.statements

def read_stmts(folder):
    fnames = get_file_names(folder)
    all_stmts = []
    for fname in fnames:
        st = get_file_stmts(fname)
        all_stmts += st
    return all_stmts

if __name__ == '__main__':
    stmts = read_stmts(base_folder)
