from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import pickle
from os import path, mkdir
from nose import SkipTest

from indra.tools.reading.read_pmids import READER_DICT, get_proc_num,\
    get_mem_total
from indra.tools.reading.read_db import _convert_id_entry, get_content,\
    get_clauses, post_reading_output, read_content, make_statements,\
    _enrich_reading_data, Reader, upload_statements

from indra.tests.test_db import get_db as get_test_db
from indra.tests.test_db import get_db_with_content

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

''' These test a depricated feature, and take FOREVER.
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


def test_reach_one_core():
    if get_mem_total() < 8:
        raise SkipTest("Not enough memory.")
    stmts = _call_reader('reach', 1)
    _check_result(stmts)


def test_reach_two_core():
    if get_mem_total() < 8:
        raise SkipTest("Not enough memory.")
    if get_proc_num() <= 2:
        raise SkipTest("Not enough processes.")
    stmts = _call_reader('reach', 2)
    _check_result(stmts)


def test_sparser_one_core():
    stmts = _call_reader('sparser', 1)
    _check_result(stmts)


def test_sparser_two_core():
    if get_proc_num() <= 2:
        raise SkipTest("Not enough processes.")
    stmts = _call_reader('sparser', 2)
    _check_result(stmts)
'''

# ==============================================================================
# Tests for NEW reading pipeline which uses the database.
# ==============================================================================


def test_convert_id_entry():
    "Test that we correctly conver the id's given us."
    id_entry = 'pmid\t: 12345\n'
    res = _convert_id_entry(id_entry)
    assert len(res) == 2 and res[0] == 'pmid' and res[1] == '12345'


def test_get_clauses():
    "Test that the clauses are correctly created."
    db = get_test_db()
    id_str_list = ['pmid:17399955']
    clauses = get_clauses(id_str_list, db.TextRef)
    assert len(clauses)


def test_get_content():
    "Test that we get content from the database successfully."
    id_str_list = ['pmid:17399955']
    content = get_content(id_str_list)
    assert len(list(content)) == 1, "Failed to get correct content."


def test_reading_content_insert():
    "Test the content primary through-put of read_db."
    db = get_db_with_content()

    print("Test reading")
    tc_list = db.select_all(db.TextContent)
    readers = [reader_class() for reader_class in Reader.__subclasses__()]
    reading_output = read_content(tc_list, readers, verbose=True,
                                  force_read=True)
    expected_output_len = len(tc_list)*len(readers)
    assert len(reading_output) == expected_output_len, \
        ("Not all text content successfully read."
         "Expected %d outputs, but got %d.") % (expected_output_len,
                                                len(reading_output))

    print("Test reading insert")
    post_reading_output(reading_output, db=db)
    r_list = db.select_all(db.Readings)

    def is_complete_match(r_list, reading_output):
        return all([any([rd.matches(r) for r in r_list])
                    for rd in reading_output.values()])

    assert is_complete_match(r_list, reading_output), \
        "Not all reading output posted."
    post_reading_output(reading_output, db=db)
    assert is_complete_match(r_list, reading_output), \
        "Uniqueness constraints failed."

    print("Test enrichement")
    assert all([rd.reading_id is None for rd in reading_output.values()]), \
        "No readings should have reading_ids already."
    _enrich_reading_data(reading_output.values(), db=db)
    assert all([rd.reading_id is not None for rd in reading_output.values()]),\
        "Some reading data objects didn't have reading_ids after enrichment."

    print("Test making statements")
    stmts = make_statements(reading_output)
    assert len(stmts), 'No statements created.'

    print("Test statement upload")
    upload_statements(stmts, db=db)
    db_stmts = db.select_all(db.Statements)
    assert len(db_stmts) == len(stmts), \
        "Only %d/%d statements added." % (len(db_stmts), len(stmts))
    assert len(db.select_all(db.Agents)), "No agents added."
