from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import pickle
from os import path, mkdir

from nose import SkipTest
from nose.plugins.attrib import attr

from indra.tools.reading.pmid_reading.read_pmids import \
    get_proc_num, get_mem_total, READER_DICT

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


def _call_reader(reader, num_cores, force_read):
    out_dir = TMP_DIR_FMT % reader
    if not path.exists(out_dir):
        mkdir(out_dir)
    stmts, pmids_unread = READER_DICT[reader](
        PMID_LIST,
        TMP_DIR_FMT % reader,
        num_cores,
        0,
        None,
        force_read,
        False
        )
    return stmts, pmids_unread


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
    assert all([stmt.evidence[0].pmid == pmid for pmid in stmts.keys()
                for stmt in stmts[pmid]]), \
        "pmids in evidence do not match true pmids."


def test_get_proc_num():
    get_proc_num()


def test_get_mem_total():
    get_mem_total()


@attr('nonpublic', 'slow')
def test_reach_one_core():
    if get_mem_total() < 8:
        raise SkipTest("Not enough memory.")
    stmts, pmids_unread = _call_reader('reach', 1, True)
    _check_result(stmts)
    assert len(pmids_unread), "Didn't read anything new."
    stmts2, pmids_unread2 = _call_reader('reach', 1, False)
    assert len(stmts) == len(stmts2)
    assert not len(pmids_unread2), "Didn't use cache."


@attr('nonpublic', 'slow')
def test_reach_two_core():
    if get_mem_total() < 8:
        raise SkipTest("Not enough memory.")
    if get_proc_num() <= 2:
        raise SkipTest("Not enough processes.")
    stmts, pmids_unread = _call_reader('reach', 2, True)
    _check_result(stmts)
    assert len(pmids_unread), "Didn't read anything new."
    stmts2, pmids_unread2 = _call_reader('reach', 2, False)
    assert len(stmts) == len(stmts2), \
        'Expected %d statements, but got %d.' % (len(stmts), len(stmts2))
    assert not len(pmids_unread2), "Didn't use cache."


@attr('nonpublic')
def test_sparser_one_core():
    stmts, pmids_unread = _call_reader('sparser', 1, True)
    _check_result(stmts)
    assert len(pmids_unread), "Didn't read anything new."
    stmts2, pmids_unread2 = _call_reader('sparser', 1, False)
    assert len(stmts) == len(stmts2)
    assert not len(pmids_unread2), "Didn't use cache."


@attr('nonpublic')
def test_sparser_two_core():
    if get_proc_num() <= 2:
        raise SkipTest("Not enough processes.")
    stmts, pmids_unread = _call_reader('sparser', 2, True)
    _check_result(stmts)
    assert len(pmids_unread), "Didn't read anything new."
    stmts2, pmids_unread2 = _call_reader('sparser', 2, False)
    assert len(stmts) == len(stmts2), \
        'Expected %d statements, but got %d.' % (len(stmts), len(stmts2))
    assert not len(pmids_unread2), "Didn't use cache."
