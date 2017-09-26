from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import json
import time
from io import BytesIO
from os import path
from docutils.io import InputError
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, UniqueConstraint, ForeignKey,\
    TIMESTAMP, create_engine, inspect, LargeBinary
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.dialects.postgresql import BYTEA
from datetime import datetime
from numbers import Number
from sqlalchemy.schema import DropTable
from sqlalchemy.ext.compiler import compiles
from indra.statements import *
from indra.util import unzip_string


@compiles(DropTable, "postgresql")
def _compile_drop_table(element, compiler, **kwargs):
    return compiler.visit_drop_table(element) + " CASCADE"


try:
    from pgcopy import CopyManager
    CAN_COPY = True
except ImportError:
    print("WARNING: pgcopy unavailable. Bulk copies will be slow.")
    CopyManager = None
    CAN_COPY = False


DEFAULTS_FILE = path.join(__path__[0], 'defaults.txt')


def _get_timestamp():
    "Get the timestamp. Needed for python 2-3 compatibility."
    try:  # Python 3
        ret = datetime.utcnow().timestamp()
    except AttributeError:  # Python 2
        now = datetime.utcnow()
        ret = time.mktime(now.timetuple())+now.microsecond/1000000.0
    return ret


def _isiterable(obj):
    "Bool determines if an object is an iterable (not a string)"
    return hasattr(obj, '__iter__') and not isinstance(obj, str)


class sqltypes:
    POSTGRESQL = 'postgresql'
    SQLITE = 'sqlite'


class texttypes:
    FULLTEXT = 'fulltext'
    ABSTRACT = 'abstract'


class formats:
    XML = 'xml'
    TEXT = 'text'
    JSON = 'json'


class DatabaseManager(object):
    """An object used to access INDRA's database.

    This object can be used to access and manage indra's database. It includes
    both basic methods and some useful, more high-level methods. It is designed
    to be used with postgresql, or sqlite.

    This object is primarily built around sqlalchemy, which is a required
    package for its use. It also optionally makes use of the pgcopy package for
    large data transfers.

    Parameters
    ----------
    host : str
        The database to which you want to interface.
    sqltype : OPTIONAL[str]
        The type of sql library used. Use one of the sql types provided by
        `sqltypes`. Default is `sqltypes.POSTGRESQL`
    """
    def __init__(self, host, sqltype=sqltypes.POSTGRESQL):
        self.host = host
        self.session = None
        self.Base = declarative_base()
        self.sqltype = sqltype

        if sqltype is sqltypes.POSTGRESQL:
            Bytea = BYTEA
        else:
            Bytea = LargeBinary

        class TextRef(self.Base):
            __tablename__ = 'text_ref'
            id = Column(Integer, primary_key=True)
            pmid = Column(String(20))
            pmcid = Column(String(20))
            doi = Column(String(100))
            pii = Column(String(250))
            url = Column(String(250), unique=True)  # Maybe longer?
            manuscript_id = Column(String(100), unique=True)
            __table_args__ = (
                UniqueConstraint('pmid', 'doi'),
                UniqueConstraint('pmcid', 'doi')
                )

        class SourceFile(self.Base):
            __tablename__ = 'source_file'
            id = Column(Integer, primary_key=True)
            source = Column(String(250), nullable=False)
            name = Column(String(250), nullable=False)
            __table_args__ = (
                UniqueConstraint('source', 'name'),
                )

        class TextContent(self.Base):
            __tablename__ = 'text_content'
            id = Column(Integer, primary_key=True)
            text_ref_id = Column(Integer,
                                 ForeignKey('text_ref.id'),
                                 nullable=False)
            text_ref = relationship(TextRef)
            source = Column(String(250), nullable=False)
            format = Column(String(250), nullable=False)
            text_type = Column(String(250), nullable=False)
            content = Column(Bytea, nullable=False)
            __table_args__ = (
                UniqueConstraint(
                    'text_ref_id', 'source', 'format', 'text_type'
                    ),
                )

        class Readings(self.Base):
            __tablename__ = 'readings'
            id = Column(Integer, primary_key=True)
            text_content_id = Column(Integer,
                                     ForeignKey('text_content.id'),
                                     nullable=False)
            text_content = relationship(TextContent)
            reader = Column(String(20), nullable=False)
            reader_version = Column(String(20), nullable=False)
            format = Column(String(20), nullable=False)  # xml, json, etc.
            bytes = Column(Bytea, nullable=False)
            ___table_args__ = (
                UniqueConstraint(
                    'text_content_id', 'reader', 'reader_version'
                    ),
                )

        class DBInfo(self.Base):
            __tablename__ = 'db_info'
            id = Column(Integer, primary_key=True)
            db_name = Column(String(20), nullable=False)
            timestamp = Column(TIMESTAMP, nullable=False)

        class Statements(self.Base):
            __tablename__ = 'statements'
            id = Column(Integer, primary_key=True)
            uuid = Column(String(20), unique=True, nullable=False)
            db_ref = Column(Integer, ForeignKey('db_info.id'))
            db_info = relationship(DBInfo)
            reader_ref = Column(Integer, ForeignKey('readings.id'))
            readings = relationship(Readings)
            type = Column(String(100), nullable=False)
            json = Column(Bytea, nullable=False)

        class Agents(self.Base):
            __tablename__ = 'agents'
            id = Column(Integer, primary_key=True)
            stmt_id = Column(Integer,
                             ForeignKey('statements.id'),
                             nullable=False)
            statements = relationship(Statements)
            db_name = Column(String(20), nullable=False)
            db_id = Column(String(20), nullable=False)
            role = Column(String(20), nullable=False)

        self.tables = {}
        for tbl in [TextRef, TextContent, Readings, SourceFile, DBInfo,
                    Statements, Agents]:
            self.tables[tbl.__tablename__] = tbl
            self.__setattr__(tbl.__name__, tbl)
        self.engine = create_engine(host)

    def __del__(self, *args, **kwargs):
        self.grab_session()
        try:
            self.session.rollback()
        except Exception as e:
            print("Failed to execute rollback of database upon deletion.")
            raise e

    def create_tables(self):
        "Create the tables for INDRA database."
        self.Base.metadata.create_all(self.engine)

    def drop_tables(self):
        "Drop all the tables for INDRA database"
        self.Base.metadata.drop_all(self.engine)

    def _clear(self):
        "Brutal clearing of all tables."
        # This is intended for testing purposes, not general use.
        # Use with care.
        self.grab_session()
        self.session.rollback()
        self.drop_tables()
        self.create_tables()

    def grab_session(self):
        "Get an active session with the database."
        if self.session is None or not self.session.is_active:
            DBSession = sessionmaker(bind=self.engine)
            self.session = DBSession()
            if self.session is None:
                raise IndraDatabaseError("Failed to grab session.")

    def get_tables(self):
        "Get a list of available tables."
        return [tbl_name for tbl_name in self.tables.keys()]

    def show_tables(self):
        "Print a list of all the available tables."
        print(self.get_tables())

    def get_active_tables(self):
        "Get the tables currently active in the database."
        return inspect(self.engine).get_table_names()

    def get_columns(self, tbl_name):
        "Get a list of the column labels for a table."
        return self.Base.metadata.tables[tbl_name].columns.keys()

    def commit(self, err_msg):
        "Commit, and give useful info if there is an exception."
        try:
            self.session.commit()
        except Exception as e:
            print(e)
            print(err_msg)
            raise

    def get_values(self, entry_list, col_names=None, keyed=False):
        "Get the column values from the entries in entry_list"
        if col_names is None and len(entry_list) > 0:  # Get everything.
            col_names = self.get_columns(entry_list[0].__tablename__)
        ret = []
        for entry in entry_list:
            if _isiterable(col_names):
                if not keyed:
                    ret.append([getattr(entry, col) for col in col_names])
                else:
                    ret.append({col: getattr(entry, col) for col in col_names})
            else:
                ret.append(getattr(entry, col_names))
        return ret

    def insert(self, tbl_name, ret_info='id', **input_dict):
        "Insert a an entry into specified table, and return id."
        self.grab_session()
        inputs = dict.fromkeys(self.get_columns(tbl_name))
        inputs.update(input_dict)
        new_entry = self.tables[tbl_name](**inputs)
        self.session.add(new_entry)
        self.commit("Excepted while trying to insert %s into %s" %
                    (inputs, tbl_name))
        return self.get_values([new_entry], ret_info)[0]

    def insert_many(self, tbl_name, input_dict_list, ret_info='id'):
        "Insert many records into the table given by table_name."
        self.grab_session()
        inputs = dict.fromkeys(self.get_columns(tbl_name))
        entry_list = []
        for input_dict in input_dict_list:
            inputs.update(input_dict)
            entry_list.append(self.tables[tbl_name](**inputs))
            inputs = inputs.fromkeys(inputs)  # Clear the values of the dict.
        self.session.add_all(entry_list)
        self.commit("Excepted while trying to insert:\n%s,\ninto %s" %
                    (input_dict_list, tbl_name))
        return self.get_values(entry_list, ret_info)

    def delete_all(self, entry_list):
        "Remove the given records from the given table."
        self.grab_session()
        for entry in entry_list:
            self.session.delete(entry)
        self.commit("Could not remove %d records from the database." %
                    len(entry_list))
        return

    def copy(self, tbl_name, data, cols=None):
        "Use pg_copy to copy over a large amount of data."
        if len(data) is 0:
            return  # Nothing to do....

        if cols is None:
            cols = self.get_columns(tbl_name)
        else:
            assert all([col in self.get_columns(tbl_name) for col in cols]),\
                "Do not recognize one of the columns in %s for table %s." % \
                (cols, tbl_name)
        if self.sqltype == sqltypes.POSTGRESQL and CAN_COPY:
            conn = self.engine.raw_connection()
            mngr = CopyManager(conn, tbl_name, cols)
            data_bts = []
            for entry in data:
                new_entry = []
                for element in entry:
                    if isinstance(element, str):
                        new_entry.append(element.encode('utf8'))
                    elif (isinstance(element, bytes)
                          or element is None
                          or isinstance(element, Number)):
                        new_entry.append(element)
                    else:
                        raise IndraDatabaseError(
                            "Don't know what to do with element of type %s."
                            "Should be str, bytes, None, or a number." %
                            type(element)
                            )
                data_bts.append(tuple(new_entry))
            mngr.copy(data_bts, BytesIO)
            conn.commit()
        else:
            # TODO: use bulk insert mappings?
            print("WARNING: You are not using postresql or do not have pgcopy,"
                  " so this will likely be very slow.")
            self.insert_many(tbl_name, [dict(zip(cols, ro)) for ro in data])

    def filter_query(self, tbls, *args):
        "Query a table and filter results."
        self.grab_session()
        if _isiterable(tbls) and not isinstance(tbls, dict):
            if isinstance(tbls[0], type(self.Base)):
                query_args = tbls
            elif isinstance(tbls[0], str):
                query_args = [self.tables[tbl] for tbl in tbls]
            else:
                raise IndraDatabaseError(
                    'Unrecognized table specification type: %s.' %
                    type(tbls[0])
                    )
        else:
            if isinstance(tbls, type(self.Base)):
                query_args = [tbls]
            elif isinstance(tbls, str):
                query_args = [self.tables[tbls]]
            else:
                raise IndraDatabaseError(
                    'Unrecognized table specification type: %s.' %
                    type(tbls)
                    )

        return self.session.query(*query_args).filter(*args)

    def select_one(self, tbls, *args):
        """Select the first value that matches requirements.

        Requirements are given in kwargs from table indicated by tbl_name. See
        *select_all*.

        Note that if your specification yields multiple results, this method
        will just return the first result without exception.
        """
        return self.filter_query(tbls, *args).first()

    def select_all(self, tbls, *args):
        """Select any and all entries from table given by tbl_name.

        The results will be filtered by your keyword arguments. For example if
        you want to get a text ref with pmid '10532205', you would call:

        .. code-block:: python

            db.select_all('text_ref', db.TextRef.pmid == '10532205')

        Note that double equals are required, not a single equal. Eqivalently
        you could call:

        .. code-block:: python

            db.select_all(db.TextRef, db.TextRef.pmid == '10532205')

        For a more complicated example, suppose you want to get all text refs
        that have full text from pmc oa, you could select:

        .. code-block:: python

           db.select_all(
               [db.TextRef, db.TextContent],
               db.TextContent.text_ref_id == db.TextRef.id,
               db.TextContent.source == 'pmc_oa',
               db.TextContent.text_type == 'fulltext'
               )
        """
        return self.filter_query(tbls, *args).all()

    def has_entry(self, tbls, *args):
        "Check whether an entry/entries matching given specs live in the db."
        q = self.filter_query(tbls, *args)
        return self.session.query(q.exists()).first()[0]

    def insert_db_stmts(self, stmts, db_name):
        "Insert statement, their database, and any affiliated agents."
        # Insert the db info
        print("Adding db %s." % db_name)
        db_ref_id = self.insert(
            'db_info',
            db_name=db_name,
            timestamp=_get_timestamp()
            )

        # Insert the statements
        for i_stmt, stmt in enumerate(stmts):
            print("Inserting stmt %s (%d/%d)" % (stmt, i_stmt+1, len(stmts)))
            stmt_id = self.insert(
                'statements',
                uuid=stmt.uuid,
                db_ref=db_ref_id,
                type=stmt.__class__.__name__,
                json=json.dumps(stmt.to_json())
                )

            # Collect the agents and add them.
            for i_ag, ag in enumerate(stmt.agent_list()):
                # If no agent, or no db_refs for the agent, skip the insert
                # that follows.
                if ag is None or ag.db_refs is None:
                    continue
                if any([isinstance(stmt, tp) for tp in
                        [Complex, SelfModification, ActiveForm]]):
                    role = 'OTHER'
                elif i_ag == 0:
                    role = 'SUBJECT'
                elif i_ag == 1:
                    role = 'OBJECT'
                else:
                    raise IndraDatabaseError("Unhandled agent role.")

                input_list = []
                for db_name, db_id in ag.db_refs.items():
                    input_list.append(
                        dict(
                            stmt_id=stmt_id,
                            role=role,
                            db_name=db_name,
                            db_id=db_id
                            )
                        )
                self.insert_many('agents', input_list)
        return

    def get_abstracts_by_pmids(self, pmid_list, unzip=True):
        "Get abstracts using the pmids in pmid_list."
        abst_list = self.filter_query(
            [self.TextRef, self.TextContent],
            self.TextContent.text_ref_id == self.TextRef.id,
            self.TextContent.text_type == 'abstract',
            self.TextRef.pmid.in_(pmid_list)
            ).all()
        if unzip:
            def unzip_func(s):
                return unzip_string(s.tobytes())
        else:
            def unzip_func(s):
                return s
        return [(r.pmid, unzip_func(c.content)) for (r, c) in abst_list]

    def get_auth_xml_pmcids(self):
        tref_list = self.filter_query(
            [self.TextRef, self.TextContent],
            self.TextRef.id == self.TextContent.text_ref_id,
            self.TextContent.text_type == texttypes.FULLTEXT,
            self.TextContent.source == 'pmc_auth'
            )
        return [tref.pmcid for tref in tref_list]

    def get_all_pmids(self):
        "Get a list of all the pmids on record."
        return self.get_values(self.select_all('text_ref'), 'pmid')

    def get_pmids(self, pmid_list):
        text_refs = self.select_all(
            'text_ref',
            self.TextRef.pmid.in_(pmid_list)
            )
        return self.get_values(text_refs, 'pmid')


class IndraDatabaseError(Exception):
    pass


def get_defaults():
    "Get the default database hosts provided in the specified `DEFAULTS_FILE`."
    default_default_file = os.path.join(
        os.path.dirname(__file__),
        'defaults.txt'
        )
    env_key_dict = {'primary': 'INDRADBPRIMARY', 'test': 'INDRADBTEST'}
    env = os.environ
    available_keys = {k: v for k, v in env_key_dict.items() if v in env.keys()}
    if not path.exists(default_default_file) and not available_keys:
        raise IndraDatabaseError(
            "Cannot find default file or environment vars."
            )
    elif path.exists(default_default_file):
        with open(default_default_file, 'r') as f:
            defaults_raw = f.read().splitlines()
        defaults_dict = {}
        for default_line in defaults_raw:
            key, value = default_line.split('=')
            defaults_dict[key.strip()] = value.strip()
    else:
        defaults_dict = {
            purpose: env_val for purpose, my_env_key in env_key_dict.items()
            for env_key, env_val in env.items() if my_env_key == env_key
            }
    return defaults_dict


def get_primary_db():
    "Get an instance to the primary database host."
    defaults = get_defaults()
    if 'primary' in defaults.keys():
        primary_host = defaults['primary']
    else:
        raise IndraDatabaseError("No primary host available in defaults file.")

    db = DatabaseManager(primary_host)
    db.grab_session()
    return db
