
from os import listdir
from indra.db import DatabaseManager

TEST_FILE = 'indra_test.db'
TEST_HOST = 'sqlite:///' + TEST_FILE

#TODO: implement setup-teardown system.

def test_create_db_file():
    "A test of the test, to ensure the file is created."
    db = DatabaseManager(TEST_HOST, sqltype='sqlite')
    db.create_tables()
    assert TEST_FILE in listdir('.'), "Test database not created."

def test_session():
    "Test whether a session can be made."
    db = DatabaseManager(TEST_HOST, sqltype='sqlite')
    db.create_tables()
    db.get_session()
    assert db.session is not None, "Session was not created."