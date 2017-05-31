from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import glob
import json
from indra import sparser
from indra.statements import stmts_from_json, get_valid_location, \
                             get_valid_residue

base_folder = os.environ['HOME'] + \
    '/data/darpa/phase3_eval/sources/sparser-20170530'

def get_file_names(base_dir):
    fnames = glob.glob(os.path.join(base_dir, '*.json'))
    return fnames

def get_file_stmts(fname):
    with open(fname, 'rb') as fh:
        print(fname)
        try:
            jd = json.load(fh)
        except ValueError as e:
            print(e)
            return []
    for st in jd:
        if st.get('type') == 'Translocation':
            for loc in ['from_location', 'to_location']:
                val = st.get(loc)
                try:
                    loc_valid = get_valid_location(val)
                    st[loc] = loc_valid
                except:
                    st[loc] = None
        try:
            res = st['residue']
            if res is False:
                st['residue'] = None
        except:
            pass

        try:
            res = st.get('residue')
            if res:
                get_valid_residue(res)
        except:
            st['residue'] = None

        try:
            res = st['position']
            if res is False:
                st['position'] = None
        except:
            pass

    stmts = stmts_from_json(jd)
    return stmts

def read_stmts(folder):
    fnames = get_file_names(folder)
    all_stmts = []
    for fname in fnames:
        st = get_file_stmts(fname)
        all_stmts += st
    return all_stmts


if __name__ == '__main__':
    stmts = read_stmts(base_folder)
