from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import pickle
from os import path, mkdir
from nose import SkipTest

from indra.tools.reading.read_db import _convert_id_entry, \
    get_content_query, get_clauses, post_reading_output, \
    get_reader_children, read_db, get_readings, _enrich_reading_data,\
    produce_statements, read_content

from indra.tests.test_db import get_db as get_test_db
from indra.tests.test_db import get_db_with_content
import random

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


def get_id_str_list(tr_list):
    id_str_list = []
    idtype_list = ['pmid', 'pmcid', 'doi']
    for tr in tr_list:
        random.shuffle(idtype_list)
        for idtype in idtype_list:
            if getattr(tr, idtype) is not None:
                break
        else:
            raise Exception("No id types found for text ref.")
        id_str_list.append('%s:%s' % (idtype, getattr(tr, idtype)))
    return id_str_list


def test_convert_id_entry():
    "Test that we correctly conver the id's given us."
    id_entry = 'pmid\t: 12345\n'
    res = _convert_id_entry(id_entry)
    assert len(res) == 2 and res[0] == 'pmid' and res[1] == '12345'


def test_get_clauses():
    "Test that the clauses are correctly created."
    db = get_test_db()
    id_str_list = ['pmid:17399955', 'pmcid:PMC3199586']
    clauses = get_clauses(id_str_list, db)
    assert len(clauses) == 1
    clause = clauses[0]
    assert 'IN' in str(clause) and 'OR' in str(clause)


def test_get_content():
    "Test that we get content from the database successfully."
    db = get_db_with_content()
    tr_list = db.select_all(db.TextRef)
    id_str_list = get_id_str_list(tr_list)
    readers = [reader_class() for reader_class in get_reader_children()]
    tc_query = get_content_query(id_str_list, readers, db=db, force_read=False)
    assert tc_query.count(),\
        "Expected some results from our query, got %d." % tc_query.count()


def test_get_reader_children():
    "Test method for getting reader objects."
    readers = get_reader_children()
    assert len(readers) == 2, \
        "Expected only 2 readers, but got %s." % str(readers)


def test_reading_content_insert():
    "Test the content primary through-put of read_db."
    db = get_db_with_content()

    print("Test reading")
    tc_list = db.select_all(db.TextContent)
    readers = [reader_class() for reader_class in get_reader_children()]
    reading_output = read_content(tc_list, readers, verbose=True)
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
                    for rd in reading_output])

    assert is_complete_match(r_list, reading_output), \
        "Not all reading output posted."
    post_reading_output(reading_output, db=db)
    assert is_complete_match(r_list, reading_output), \
        "Uniqueness constraints failed."

    print("Test enrichement")
    assert all([rd.reading_id is None for rd in reading_output]), \
        "No readings should have reading_ids already."
    _enrich_reading_data(reading_output, db=db)
    assert all([rd.reading_id is not None for rd in reading_output]),\
        "Some reading data objects didn't have reading_ids after enrichment."

    print("Test making statements")
    stmts = produce_statements(reading_output, db=db)
    assert len(stmts), 'No statements created.'
    db_stmts = db.select_all(db.Statements)
    assert len(db_stmts) == len(stmts), \
        "Only %d/%d statements added." % (len(db_stmts), len(stmts))
    assert len(db.select_all(db.Agents)), "No agents added."


def test_read_db():
    "Test the read_db function with various settings."
    # Prep the inputs.
    db = get_db_with_content()
    complete_tr_list = db.select_all(db.TextRef)
    id_str_list = get_id_str_list(complete_tr_list)
    readers = [reader_class() for reader_class in get_reader_children()
               if reader_class.name == 'SPARSER']

    # Run the reading with default batch size, no force_fulltext, and
    # no force_read
    reading_output_1 = read_db(id_str_list, readers, db=db)
    N1 = len(reading_output_1)
    N1_exp = len(readers)*db.filter_query(db.TextContent).count()
    assert N1 == N1_exp, \
        'Expected %d readings, but got %d.' % (N1_exp, N1)
    post_reading_output(reading_output_1, db=db)  # setup for later test.
    N1_db = len(db.select_all(db.Readings))
    assert N1_db == N1, \
        'Expected %d readings to be copied to db, only %d found.' % (N1, N1_db)

    # Run the reading with default batch size, no force_fulltext, but with
    # force_read = True (this hould produce new readings.)
    reading_output_2 = read_db(id_str_list, readers, db=db,
                               force_read=True)
    N2 = len(reading_output_2)
    assert N1 == N2, "Got %d readings from run 1 but %d from run 2." % (N1, N2)

    # Run the reading with default batch size, no force_fulltext, but without
    # force_read = True (this should NOT produce new readings.)
    old_readings = get_readings(id_str_list, readers, db=db)
    reading_output = read_db(id_str_list, readers, db=db)
    assert len(reading_output) == 0, "Got new readings when force_read=False."
    assert len(old_readings) == N1, \
        "Did not get old readings when force_read=False."
