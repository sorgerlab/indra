from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import pickle
import random
import zlib
from os import path
from nose.plugins.attrib import attr

from indra.tools.reading.db_reading import read_db as rdb
from indra.tools.reading.read_files import read_files
from indra.tools.reading.util.script_tools import make_statements
from indra.tools.reading.readers import SparserReader
from indra.tools.reading.readers import get_readers as get_all_readers

from indra.db import formats
from indra.tests.test_db import get_db as get_test_db
from indra.tests.test_db import get_db_with_pubmed_content


# ==============================================================================
# Tests for NEW reading pipeline which uses the database.
# ==============================================================================


def get_id_dict(tr_list):
    idtype_list = ['pmid', 'pmcid', 'doi']
    id_dict = {id_type: [] for id_type in idtype_list}
    for tr in tr_list:
        random.shuffle(idtype_list)
        for idtype in idtype_list:
            if getattr(tr, idtype) is not None:
                break
        else:
            raise Exception("No id types found for text ref.")
        id_dict[idtype].append(getattr(tr, idtype))
    return id_dict


def get_readers(*names, **kwargs):
    return [reader_class(**kwargs) for reader_class in get_all_readers()
            if (not names or reader_class.name in names)]


def test_convert_id_entry():
    "Test that we correctly conver the id's given us."
    id_entry = 'pmid\t: 12345\n'
    res = rdb._convert_id_entry(id_entry)
    assert len(res) == 2 and res[0] == 'pmid' and res[1] == '12345'


def test_get_clauses():
    "Test that the clauses are correctly created."
    db = get_test_db()
    id_dict = {'pmid': '17399955', 'pmcid': 'PMC3199586'}
    clauses = rdb.get_clauses(id_dict, db)
    assert len(clauses) == 2
    clause = clauses[0]
    assert 'IN' in str(clause),\
        "Unexpected form for clause: %s" % str(clause)


@attr('nonpublic')
def test_get_content():
    "Test that content querries are correctly formed."
    db = get_db_with_pubmed_content()
    tr_list = db.select_all(db.TextRef)
    id_dict = get_id_dict(tr_list)
    readers = get_readers()
    tc_query_1 = rdb.get_content_query(id_dict, readers, db=db,
                                       force_read=False)
    N_exp = db.filter_query(db.TextContent).count()
    N_1 = tc_query_1.count()
    assert N_1 == N_exp,\
        "Expected %d results in our query, got %d." % (N_exp, N_1)

    # This tests both that tcid and trid are recognized, and that we really
    # are getting just a subset of the content.
    test_tcs_2 = tc_query_1.limit(2).all()
    small_id_dict = {'tcid': [test_tcs_2[0].id],
                     'trid': [test_tcs_2[1].text_ref_id]}
    tc_query_2 = rdb.get_content_query(small_id_dict, readers, db=db)
    N_2 = tc_query_2.count()
    assert N_2 == 2, "Expected 2 items in query, but got %d." % N_2

    # Now test the 'all' feature.
    tc_query_3 = rdb.get_content_query('all', readers, db=db)
    N_3 = tc_query_3.count()
    assert N_3 == N_1, \
        "Expected to get %d items in query, but got %d." % (N_1, N_3)

    # Test response to empyt dict.
    assert rdb.get_content_query({}, readers, db=db) is None, \
        "Expected None when passing no ids."


@attr('nonpublic')
def test_get_reader_children():
    "Test method for getting reader objects."
    readers = get_readers()
    assert len(readers) == 2, \
        "Expected only 2 readers, but got %s." % str(readers)


@attr('slow', 'nonpublic')
def test_reading_content_insert():
    "Test the content primary through-put of make_db_readings."
    db = get_db_with_pubmed_content()

    print("Test reading")
    tc_list = db.select_all(db.TextContent)
    readers = get_readers()
    reading_output = []
    for reader in readers:
        reading_output += reader.read(tc_list, verbose=True)
    expected_output_len = len(tc_list)*len(readers)
    assert len(reading_output) == expected_output_len, \
        ("Not all text content successfully read."
         "Expected %d outputs, but got %d.") % (expected_output_len,
                                                len(reading_output))

    print("Test reading insert")
    rdb.upload_readings(reading_output, db=db)
    r_list = db.select_all(db.Readings)

    def is_complete_match(r_list, reading_output):
        return all([any([rd.matches(r) for r in r_list])
                    for rd in reading_output])

    assert is_complete_match(r_list, reading_output), \
        "Not all reading output posted."
    rdb.upload_readings(reading_output, db=db)
    assert is_complete_match(r_list, reading_output), \
        "Uniqueness constraints failed."

    print("Test enrichement")
    assert all([rd.reading_id is None for rd in reading_output]), \
        "No readings should have reading_ids already."
    rdb._enrich_reading_data(reading_output, db=db)
    assert all([rd.reading_id is not None for rd in reading_output]),\
        "Some reading data objects didn't have reading_ids after enrichment."

    print("Test making statements")
    stmts = rdb.produce_statements(reading_output, db=db)
    assert len(stmts), 'No statements created.'
    db_stmts = db.select_all(db.Statements)
    assert len(db_stmts) == len(stmts), \
        "Only %d/%d statements added." % (len(db_stmts), len(stmts))
    assert len(db.select_all(db.Agents)), "No agents added."


@attr('nonpublic')
def test_read_db():
    "Test the low level make_db_readings functionality with various settings."
    # Prep the inputs.
    db = get_db_with_pubmed_content()
    complete_tr_list = db.select_all(db.TextRef)
    id_dict = get_id_dict(complete_tr_list)
    readers = get_readers('SPARSER')

    # Run the reading with default batch size, no force_fulltext, and
    # no force_read
    reading_output_1 = rdb.make_db_readings(id_dict, readers, db=db)
    N1 = len(reading_output_1)
    N1_exp = len(readers)*db.filter_query(db.TextContent).count()
    assert N1 == N1_exp, \
        'Expected %d readings, but got %d.' % (N1_exp, N1)
    rdb.upload_readings(reading_output_1, db=db)  # setup for later test.
    N1_db = len(db.select_all(db.Readings))
    assert N1_db == N1, \
        'Expected %d readings to be copied to db, only %d found.' % (N1, N1_db)

    # Run the reading with default batch size, no force_fulltext, but with
    # force_read = True (this hould produce new readings.)
    reading_output_2 = rdb.make_db_readings(id_dict, readers, db=db,
                                            force_read=True)
    N2 = len(reading_output_2)
    assert N1 == N2, "Got %d readings from run 1 but %d from run 2." % (N1, N2)

    # Run the reading with default batch size, no force_fulltext, but without
    # force_read = True (this should NOT produce new readings.)
    old_readings = rdb.get_db_readings(id_dict, readers, db=db)
    reading_output = rdb.make_db_readings(id_dict, readers, db=db)
    assert len(reading_output) == 0, "Got new readings when force_read=False."
    assert len(old_readings) == N1, \
        "Did not get old readings when force_read=False."


@attr('slow', 'nonpublic')
def test_produce_readings():
    "Comprehensive test of the high level production of readings."
    # Prep the inputs.
    db = get_db_with_pubmed_content()
    complete_tr_list = db.select_all(db.TextRef)
    id_dict = get_id_dict(complete_tr_list)

    # Test with just sparser for tollerable speeds.
    reader_list = get_readers('SPARSER')

    # Test the read_mode='none' option (should yield nothing, because there
    # aren't any readings yet.)
    outputs_0 = rdb.produce_readings(id_dict, reader_list, verbose=True, db=db,
                                     read_mode='none')
    assert len(outputs_0) == 0

    # Test just getting a pickle file (Nothing should be posted to the db.).
    pkl_file = 'test_db_res.pkl'
    outputs_1 = rdb.produce_readings(id_dict, reader_list, verbose=True, db=db,
                                     no_upload=True, pickle_file=pkl_file)
    N_out = len(outputs_1)
    N_exp = len(reader_list)*db.filter_query(db.TextContent).count()
    assert N_out == N_exp, "Expected %d readings, got %d." % (N_exp, N_out)
    assert path.exists(pkl_file), "Pickle file not created."
    with open(pkl_file, 'rb') as f:
        N_pkl = len(pickle.load(f))
    assert N_pkl == N_exp, \
        "Expected %d readings in pickle, got %d." % (N_exp, N_out)
    N_readings = db.filter_query(db.Readings).count()
    assert N_readings == 0, \
        "There shouldn't be any readings yet, but found %d." % N_readings

    # Test reading and insert to the database.
    rdb.produce_readings(id_dict, reader_list, verbose=True, db=db)
    N_db = db.filter_query(db.Readings).count()
    assert N_db == N_exp, "Excpected %d readings, got %d." % (N_exp, N_db)

    # Test reading again, without read_mode='all'
    outputs_2 = rdb.produce_readings(id_dict, reader_list, verbose=True, db=db)
    assert len(outputs_2) == N_exp, \
        "Got %d readings, expected %d." % (len(outputs_2), N_exp)
    assert all([rd.reading_id is not None for rd in outputs_2])

    # Test with read_mode='none' again.
    outputs_3 = rdb.produce_readings(id_dict, reader_list, verbose=True, db=db,
                                     read_mode='none')
    assert len(outputs_3) == N_exp
    assert all([rd.reading_id is not None for rd in outputs_3])

    # Test the read_mode='all'.
    outputs_4 = rdb.produce_readings(id_dict, reader_list, verbose=True, db=db,
                                     read_mode='all')
    assert len(outputs_4) == N_exp
    assert all([rd.reading_id is None for rd in outputs_4])


@attr('slow', 'nonpublic')
def test_read_files():
    "Test that the system can read files."
    db = get_db_with_pubmed_content()

    # Create the test files.
    test_file_fmt = 'test_reading_input.%s'
    example_files = []
    for fmt in [formats.TEXT, formats.XML]:
        tc = db.select_one(db.TextContent, db.TextContent.format == fmt)
        if tc is None:
            print("Could not find %s content for testing." % fmt)
            continue
        suffix = fmt
        if fmt is formats.XML:
            suffix = 'n' + fmt
        with open(test_file_fmt % suffix, 'wb') as f:
            f.write(zlib.decompress(tc.content, 16+zlib.MAX_WBITS))
        example_files.append(test_file_fmt % suffix)

    # Now read them.
    readers = get_readers()
    outputs = read_files(example_files, readers)
    N_out = len(outputs)
    N_exp = len(example_files)
    assert N_out == N_exp, "Expected %d outputs, got %d." % (N_exp, N_out)


@attr('nonpublic')
def test_sparser_parallel():
    "Test running sparser in parallel."
    db = get_db_with_pubmed_content()
    sparser_reader = SparserReader(n_proc=2)
    tc_list = db.select_all(db.TextContent)
    result = sparser_reader.read(tc_list, verbose=True, log=True)
    N_exp = len(tc_list)
    N_res = len(result)
    assert N_exp == N_res, \
        "Expected to get %d results, but got %d." % (N_exp, N_res)


@attr('nonpublic')
def test_sparser_parallel_one_batch():
    "Test that sparser runs with multiple procs with batches of 1."
    db = get_db_with_pubmed_content()
    sparser_reader = SparserReader(n_proc=2)
    tc_list = db.select_all(db.TextContent)
    result = sparser_reader.read(tc_list, verbose=True, n_per_proc=1)
    N_exp = len(tc_list)
    N_res = len(result)
    assert N_exp == N_res, \
        "Expected to get %d results, but got %d." % (N_exp, N_res)


@attr('slow', 'nonpublic')
def test_multi_batch_run():
    "Test that reading works properly with multiple batches run."
    db = get_db_with_pubmed_content()
    readers = get_readers()
    tc_list = db.select_all(db.TextContent)
    id_dict = {'tcid': [tc.id for tc in tc_list]}
    outputs = rdb.make_db_readings(id_dict, readers,
                                   batch_size=len(tc_list)//2, db=db)
    # This should catch any repeated readings.
    rdb.upload_readings(outputs, db=db)
    num_readings = db.filter_query(db.Readings).count()
    assert num_readings == 2*len(tc_list), \
        "Expected %d readings, only found %d." % (2*len(tc_list), num_readings)


@attr('slow', 'nonpublic')
def test_multiproc_statements():
    "Test the multiprocessing creation of statements."
    db = get_db_with_pubmed_content()
    readers = get_readers()
    tc_list = db.select_all(db.TextContent)
    id_dict = {'tcid': [tc.id for tc in tc_list]}
    outputs = rdb.make_db_readings(id_dict, readers, db=db)
    stmts = make_statements(outputs, 2)
    assert len(stmts)
