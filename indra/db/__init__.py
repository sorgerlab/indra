from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['sqltypes', 'texttypes', 'formats', 'DatabaseManager',
           'IndraDatabaseError', 'get_defaults', 'get_primary_db']

import json
import time
import re
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
from sqlalchemy.sql.expression import Delete, Update
from sqlalchemy.ext.compiler import compiles
from indra.statements import *
from indra.util import unzip_string
from indra.util.get_version import get_version


# Solution to fix postgres drop tables
# See: https://stackoverflow.com/questions/38678336/sqlalchemy-how-to-implement-drop-table-cascade
@compiles(DropTable, "postgresql")
def _compile_drop_table(element, compiler, **kwargs):
    return compiler.visit_drop_table(element) + " CASCADE"


# Solution to fix deletes with constraints from multiple tables.
# See: https://groups.google.com/forum/#!topic/sqlalchemy/cIvgH2y01_o
@compiles(Delete)
def compile_delete(element, compiler, **kw):
    text = compiler.visit_delete(element, **kw)
    extra_froms = Update._extra_froms.__get__(element)
    if extra_froms:
        text = re.sub(
                    r"(FROM \S+)",
                    lambda m: "%s USING %s" % (
                        m.group(1),
                        ", ".join(
                            compiler.process(fr, asfrom=True, **kw)
                            for fr in extra_froms
                        )
                    ),
                    text
                )
    return text


try:
    from pgcopy import CopyManager
    CAN_COPY = True
except ImportError:
    print("WARNING: pgcopy unavailable. Bulk copies will be slow.")
    CopyManager = None
    CAN_COPY = False


DEFAULTS_FILE = path.join(__path__[0], 'defaults.txt')
__PRIMARY_DB = None


def _isiterable(obj):
    "Bool determines if an object is an iterable (not a string)"
    return hasattr(obj, '__iter__') and not isinstance(obj, str)


class _map_class(object):
    @classmethod
    def _getattrs(self):
        return {
            k: v for k, v in self.__dict__.items() if not k.startswith('_')
            }

    @classmethod
    def items(self):
        return self._getattrs().items()

    @classmethod
    def values(self):
        return self._getattrs().values()

    @classmethod
    def keys(self):
        return self._getattrs().keys()


class sqltypes(_map_class):
    POSTGRESQL = 'postgresql'
    SQLITE = 'sqlite'


class texttypes(_map_class):
    FULLTEXT = 'fulltext'
    ABSTRACT = 'abstract'


class formats(_map_class):
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

    If you wish to access the primary database, you can simply use the
    `get_primary_db` to get an instance of this object using the default
    settings.

    Parameters
    ----------
    host : str
        The database to which you want to interface.
    sqltype : OPTIONAL[str]
        The type of sql library used. Use one of the sql types provided by
        `sqltypes`. Default is `sqltypes.POSTGRESQL`
    label : OPTIONAL[str]
        A short string to indicate the purpose of the db instance. Set as
        primary when initialized be `get_primary_db`.

    Example
    -------
    If you wish to acces the primary database and find the the metadata for a
    particular pmid, 1234567:

    >> from indra.db import get_primary_db()
    >> db = get_primary_db()
    >> res = db.select_all(db.TextRef, db.TextRef.pmid == '1234567')

    You will get a list of objects whose attributes give the metadata contained
    in the columns of the table.

    For more sophisticated examples, several use cases can be found in
    `indra.tests.test_db`.
    """
    def __init__(self, host, sqltype=sqltypes.POSTGRESQL, label=None):
        self.host = host
        self.session = None
        self.Base = declarative_base()
        self.sqltype = sqltype
        self.label = label

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
            uuid = Column(String(40), unique=True, nullable=False)
            db_ref = Column(Integer, ForeignKey('db_info.id'))
            db_info = relationship(DBInfo)
            reader_ref = Column(Integer, ForeignKey('readings.id'))
            readings = relationship(Readings)
            type = Column(String(100), nullable=False)
            indra_version = Column(String(100), nullable=False)
            json = Column(Bytea, nullable=False)

        class Agents(self.Base):
            __tablename__ = 'agents'
            id = Column(Integer, primary_key=True)
            stmt_id = Column(Integer,
                             ForeignKey('statements.id'),
                             nullable=False)
            statements = relationship(Statements)
            db_name = Column(String(40), nullable=False)
            db_id = Column(String(40), nullable=False)
            role = Column(String(20), nullable=False)

        self.tables = {}
        for tbl in [TextRef, TextContent, Readings, SourceFile, DBInfo,
                    Statements, Agents]:
            self.tables[tbl.__tablename__] = tbl
            self.__setattr__(tbl.__name__, tbl)
        self.engine = create_engine(host)

    def __del__(self, *args, **kwargs):
        try:
            self.grab_session()
            self.session.rollback()
        except Exception as e:
            print("Failed to execute rollback of database upon deletion.")
            raise e

    def create_tables(self, tbl_list=None):
        "Create the tables for INDRA database."
        if tbl_list is None:
            self.Base.metadata.create_all(self.engine)
        else:
            for tbl in tbl_list:
                if isinstance(tbl, str):
                    tbl = self.tables[tbl]
                tbl.__table__.create(bind=self.engine)
        return

    def drop_tables(self, tbl_list=None, force=False):
        """Drop the tables for INDRA database given in tbl_list.

        If tbl_list is None, all tables will be dropped. Note that if `force`
        is False, a warning prompt will be raised to asking for confirmation,
        as this action will remove all data from that table.
        """
        if not force:
            if tbl_list is None:
                msg = "Do you really want to clear the primary database? [y/N]: "
                try:
                    resp = raw_input(msg)
                except NameError:
                    resp = input(msg)
                if resp != 'y' and resp != 'yes':
                    logger.info('Aborting clearing of database.')
                    return
            else:
                msg = "You are going to clear the following tables:\n"
                msg += str(tbl_list) + '\n'
                msg += "Do you really want to clear these tables? [y/N]: "
        if tbl_list is None:
            self.Base.metadata.drop_all(self.engine)
        else:
            for tbl in tbl_list:
                if isinstance(tbl, str):
                    tbl = self.tables[tbl]
                tbl.__table__.drop()
        return

    def _clear(self, tbl_list=None, force=False):
        "Brutal clearing of all tables in tbl_list, or all tables."
        # This is intended for testing purposes, not general use.
        # Use with care.
        self.grab_session()
        self.session.rollback()
        self.drop_tables(tbl_list, force=force)
        self.create_tables(tbl_list)
        return

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
        if isinstance(tbl_name, type(self.Base)):
            tbl_name = tbl_name.__tablename__
        return self.Base.metadata.tables[tbl_name].columns.keys()

    def commit(self, err_msg):
        "Commit, and give useful info if there is an exception."
        try:
            self.session.commit()
        except Exception as e:
            if self.session is not None:
                self.session.rollback()
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
            logger.warning("You are not using postgresql or do not have "
                           "pgcopy, so this will likely be very slow.")
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

    def insert_db_stmts(self, stmts, db_ref_id):
        "Insert statement, their database, and any affiliated agents."
        # Prepare the statements for copying
        stmt_data = []
        cols = ('uuid', 'db_ref', 'type', 'json', 'indra_version')
        for stmt in stmts:
            stmt_rec = (
                stmt.uuid,
                db_ref_id,
                stmt.__class__.__name__,
                json.dumps(stmt.to_json()).encode('utf8'),
                get_version()
            )
            stmt_data.append(stmt_rec)
        self.copy('statements', stmt_data, cols)

        # Build a dict mapping stmt UUIDs to statement IDs
        uuid_list = [s.uuid for s in stmts]
        stmt_rec_list = self.select_all('statements',
                                        self.Statements.uuid.in_(uuid_list))
        stmt_uuid_dict = {uuid: sid for uuid, sid in
                          self.get_values(stmt_rec_list, ['uuid', 'id'])}

        # Now assemble agent records
        agent_data = []
        for stmt in stmts:
            stmt_id = stmt_uuid_dict[stmt.uuid]
            for ag_ix, ag in enumerate(stmt.agent_list()):
                # If no agent, or no db_refs for the agent, skip the insert
                # that follows.
                if ag is None or ag.db_refs is None:
                    continue
                if any([isinstance(stmt, tp) for tp in
                        [Complex, SelfModification, ActiveForm]]):
                    role = 'OTHER'
                elif ag_ix == 0:
                    role = 'SUBJECT'
                elif ag_ix == 1:
                    role = 'OBJECT'
                else:
                    raise IndraDatabaseError("Unhandled agent role.")
                for ns, ag_id in ag.db_refs.items():
                    ag_rec = (stmt_id, ns, ag_id, role)
                    agent_data.append(ag_rec)
        cols = ('stmt_id', 'db_name', 'db_id', 'role')
        self.copy('agents', agent_data, cols)
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

    def _get_timestamp(self):
        "Get the timestamp. Needed for python 2-3 compatibility."
        try:  # Python 3
            ret = datetime.utcnow()
            #if self.sqltype != sqltypes.SQLITE:
            #    ret = ret.timestamp()
        except AttributeError:  # Python 2
            ret = datetime.utcnow()
            #if self.sqltype != sqltypes.SQLITE:
            #    ret = time.mktime(ret.timetuple())+ret.microsecond/1000000.0
        return ret


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


def get_primary_db(force_new=False):
    """Get a DatabaseManager instance for the primary database host.

    The primary database host is defined in the defaults.txt file, or in a file
    given by the environment variable DEFAULTS_FILE. Alternatively, it may be
    defined by the INDRADBPRIMARY environment variable. If none of the above
    are specified, this function will raise an exception.

    Parameters
    ----------
    force_new : bool
        If true, a new instance will be created and returned, regardless of
        whether there is an existing instance or not. Default is False, so that
        if this function has been called before within the global scope, a the
        instance that was first created will be returned.

    Returns
    -------
    primary_db : DatabaseManager instance
        An instance of the database manager that is attached to the primary
        database.

    Note: by default, calling this function twice will return the same
    `DatabaseManager` instance. In other words:

    > db1 = get_primary_db()
    > db2 = get_primary_db()
    > db1 is db2
    True

    This means also that, for example `db1.select_one(db2.TextRef)` will work,
    in the above context.

    It is still recommended that when creating a script or function, or other
    general application, you should not rely on this feature to get your access
    to the database, as it can make substituting a different database host both
    complicated and messy. Rather, a database instance should be explicitly
    passed between different users as is done in the `by_gene_role_type`
    function's call to `get_statements` in `indra.db.query_db_stmts`.
    """
    defaults = get_defaults()
    if 'primary' in defaults.keys():
        primary_host = defaults['primary']
    else:
        raise IndraDatabaseError("No primary host available in defaults file.")

    global __PRIMARY_DB
    if __PRIMARY_DB is None or force_new:
        __PRIMARY_DB = DatabaseManager(primary_host, label='primary')
        __PRIMARY_DB.grab_session()
    return __PRIMARY_DB
