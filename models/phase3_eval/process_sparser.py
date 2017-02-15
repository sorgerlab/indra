from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import glob
from indra import sparser

base_folder = os.path.join(os.environ['HOME'],
                           'data/darpa/phase3_eval/sources/sparser-20170210')

def get_file_names(base_dir):
    fnames = glob.glob(os.path.join(base_dir, '*.xml'))
    return fnames

def get_file_stmts(fname):
    with open(fname, 'rb') as fh:
        xml_bytes = fh.read()
        xml_bytes = xml_bytes.replace(b'hmsid', b'pmid')
        sp = sparser.process_xml(xml_bytes)
        if sp is None:
            print('ERROR: Could not process %s' % fname.split('/')[-1])
            print('----')
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
