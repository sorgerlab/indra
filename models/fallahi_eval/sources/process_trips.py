from indra.util import _require_python3
import os
import glob
import pickle
from indra.sources import trips
from util import prefixed_pkl, based


base_folder = os.path.join(based, 'trips-20170206')


def get_file_names(base_dir):
    fnames = glob.glob(os.path.join(base_dir, '*.ekb'))
    return fnames


def get_file_stmts(fname):
    with open(fname, 'rt') as fh:
        xml_str = fh.read()
        tp = trips.process_xml(xml_str)
        if tp is None:
            return []
        return tp.statements


def read_stmts(folder):
    fnames = get_file_names(folder)
    all_stmts = []
    for i, fname in enumerate(fnames):
        print('%d/%d' % (i, len(fnames)))
        print(fname)
        print('='*len(fname))
        st = get_file_stmts(fname)
        all_stmts += st
    return all_stmts


if __name__ == '__main__':
    stmts = read_stmts(base_folder)
    print('Collected %d Statements from TRIPS' % len(stmts))
    with open(prefixed_pkl('trips'), 'wb') as fh:
        pickle.dump(stmts, fh)
