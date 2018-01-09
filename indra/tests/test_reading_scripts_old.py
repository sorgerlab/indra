from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import pickle
from os import path, mkdir

from nose import SkipTest
from nose.plugins.attrib import attr

from indra.tools.reading.read_pmids import get_proc_num, get_mem_total,\
    READER_DICT

# ==============================================================================
# Tests for OLD reading pipeline that did not use the database.
# ==============================================================================


PMID_LIST = [
    '27085964',
    '28739733',
    '18201725',
    '21655183',
    '16254254',
    '15899863',
    '12122017',
    '19426868',
    # '18317068',
    # '24657168'
    ]
BASENAME = 'test_tmp'
TMP_DIR_FMT = '%s_%%s' % BASENAME
OUTPUT_FILE_FMT = '%s_stmts_0-10.pkl' % TMP_DIR_FMT
READINGS_PKL = 'sample_reach_outputs.pkl'


def _call_reader(reader, num_cores):
    out_dir = TMP_DIR_FMT % reader
    if not path.exists(out_dir):
        mkdir(out_dir)
    stmts = READER_DICT[reader](
        PMID_LIST,
        TMP_DIR_FMT % reader,
        num_cores,
        0,
        len(PMID_LIST),
        True,
        False
        )
    return stmts


def _check_blind_result(reader):
    output_file = OUTPUT_FILE_FMT % reader
    assert path.exists(output_file),\
        "Expected output pickle file missing: %s." % output_file
    with open(output_file, 'rb') as f:
        pkl_out = pickle.load(f)
    assert reader in pkl_out.keys(),\
        "Pickle file does not contain key for reader."
    assert len(pkl_out[reader]),\
        "No statements found."


def _check_result(stmts):
    assert len(stmts), "No statements found."


def test_get_proc_num():
    get_proc_num()


def test_get_mem_total():
    get_mem_total()


@attr('nonpublic', 'slow')
def test_reach_one_core():
    if get_mem_total() < 8:
        raise SkipTest("Not enough memory.")
    stmts = _call_reader('reach', 1)
    _check_result(stmts)


@attr('nonpublic', 'slow')
def test_reach_two_core():
    if get_mem_total() < 8:
        raise SkipTest("Not enough memory.")
    if get_proc_num() <= 2:
        raise SkipTest("Not enough processes.")
    stmts = _call_reader('reach', 2)
    _check_result(stmts)


@attr('nonpublic')
def test_sparser_one_core():
    stmts = _call_reader('sparser', 1)
    _check_result(stmts)


@attr('nonpublic')
def test_sparser_two_core():
    if get_proc_num() <= 2:
        raise SkipTest("Not enough processes.")
    stmts = _call_reader('sparser', 2)
    _check_result(stmts)
