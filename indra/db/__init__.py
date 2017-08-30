from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import json
from indra.statements import *
from indra.util import unzip_string
from indra.databases import hgnc_client
from io import BytesIO
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, UniqueConstraint, ForeignKey,\
    TIMESTAMP, create_engine, inspect
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.dialects.postgresql import BYTEA
from datetime import datetime
import time
try:
    from pgcopy import CopyManager
    CAN_COPY = True
except:
    print("WARNING: pgcopy unavailable. Bulk copies will be slow.")
    CopyManager = None
    CAN_COPY = False

from docutils.io import InputError
#from lxml.includes.xpath import XPATH_INVALID_TYPE

DEFAULT_AWS_HOST = 'indradb.cwcetxbvbgrf.us-east-1.rds.amazonaws.com'


def get_timestamp():
    "Get the timestamp. Needed for python 2-3 compatibility."
    try: # Python 3
        ret = datetime.utcnow().timestamp()
    except: # Python 2
        now = datetime.utcnow()
        ret = time.mktime(now.timetuple())+now.microsecond/1000000.0
    return ret


def isiterable(obj):
    "Bool determines if an object is an iterable (not a string)"
    return hasattr(obj, '__iter__') and not isinstance(obj, str)


class sqltypes:
    POSTGRESQL = 'postgresql'
    SQLITE = 'sqlite'


class texttypes:
    FULLTEXT = 'fulltext'
    ABSTRACT = 'abstract'


class DatabaseManager(object):
    """An object used to access indra's database.
    
    This object can be used to access and manage indra's database. It includes
    both basic methods and some useful, more high-level methods. It is designed
    to be used with postresql, however for testing, it uses a local sqlite 
    database.
    
    This object is primarily built around sqlalchemy, which is a required
    package for its use. It also optionally makes use of the pgcopy package for
    large data transfers.
    
    Init:
    ----
    host -- the database to which you want to interface.
    sqltype (optional) -- the type of sql library used. Use one of the sql
            types provided by `sqltypes`.
    """
    def __init__(self, host, sqltype=sqltypes.POSTGRESQL):
        self.host = host
        self.Base = declarative_base()
        self.sqltype = sqltype
        
        if sqltype is sqltypes.POSTGRESQL:
            Bytea = BYTEA
        else:
            Bytea = String
        
        class TextRef(self.Base):
            __tablename__ = 'text_ref'
            id = Column(Integer, primary_key=True)
            pmid = Column(String(20))
            pmcid = Column(String(20))
            doi = Column(String(100))
            pii = Column(String(250))
            url = Column(String(250), unique = True) # Maybe longer?
            manuscript_id = Column(String(100), unique = True)
            __table_args__ = (
                UniqueConstraint('pmid', 'doi'),
                UniqueConstraint('pmcid', 'doi')
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
                UniqueConstraint('text_ref_id', 'source', 'format', 'text_type'),
                )
        
        class Reach(self.Base):
            __tablename__ = 'reach'
            id = Column(Integer, primary_key=True)
            text_content_id = Column(Integer, 
                                     ForeignKey('text_content.id'), 
                                     nullable=False)
            text_content = relationship(TextContent)
            version = Column(String(20), nullable=False)
            json = Column(Bytea, nullable=False)
            ___table_args__ = (
                UniqueConstraint('text_content_id', 'version'),
                )
        
        class DBInfo(self.Base):
            __tablename__ = 'db_info'
            id = Column(Integer, primary_key=True)
            db_name = Column(String(20), nullable=False)
            #FIXME: I don't think the default here is really correct.
            timestamp = Column(TIMESTAMP, nullable=False)
        
        class Statements(self.Base):
            __tablename__ = 'statements'
            id = Column(Integer, primary_key=True)
            uuid = Column(String(20), unique=True, nullable=False)
            db_ref = Column(Integer, ForeignKey('db_info.id'))
            db_info = relationship(DBInfo)
            type = Column(String(100), nullable=False)
            json = Column(String(500), nullable=False)
        
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
        for tbl in [TextRef, TextContent, Reach, DBInfo, Statements, Agents]:
            self.tables[tbl.__tablename__] = tbl
            self.__setattr__(tbl.__name__, tbl)
        self.engine = create_engine(host)
        self.session = None


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
        self.drop_tables()
        self.create_tables()


    def grab_session(self):
        "Get an active session with the database."
        if self.session is None or not self.session.is_active:
            DBSession = sessionmaker(bind=self.engine)
            self.session = DBSession()


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
        except:
            print(err_msg)
            raise


    def get_values(self, entry_list, col_names=None, keyed=False):
        "Get the column values from the entries in entry_list"
        if col_names is None and len(entry_list) > 0: # Get everything.
            col_names = self.get_columns(entry_list[0].__tablename__)
        ret = []
        for entry in entry_list:
            if isiterable(col_names):
                if not keyed:
                    ret.append([getattr(entry, col) for col in col_names])
                else:
                    ret.append({col:getattr(entry, col) for col in col_names})
            else:
                ret.append(getattr(entry, col_names))
        return ret 


    def insert(self, tbl_name, ret_info = 'id', **input_dict):
        "Insert a an entry into specified table, and return id."
        self.grab_session()
        inputs = dict.fromkeys(self.get_columns(tbl_name))
        inputs.update(input_dict)
        new_entry = self.tables[tbl_name](**inputs)
        self.session.add(new_entry)
        self.commit("Excepted while trying to insert %s into %s" % 
                  (inputs, tbl_name))
        return self.get_values([new_entry], ret_info)[0]


    def insert_many(self, tbl_name, input_dict_list, ret_info = 'id'):
        "Insert many records into the table given by table_name."
        self.grab_session()
        inputs = dict.fromkeys(self.get_columns(tbl_name))
        entry_list = []
        for input_dict in input_dict_list:
            inputs.update(input_dict)
            entry_list.append(self.tables[tbl_name](**inputs))
            inputs = inputs.fromkeys(inputs) # Clear the values of the dict.
        self.session.add_all(entry_list)
        self.commit("Excepted while trying to insert:\n%s,\ninto %s" % 
                  (input_dict_list, tbl_name))
        return self.get_values(entry_list, ret_info)


    def copy(self, tbl_name, data, cols=None):
        "Use pg_copy to copy over a large amount of data."
        if cols is None:
            cols = self.get_columns(tbl_name)
        else:
            assert all([col in self.get_columns(tbl_name) for col in cols]),\
                "Do not recognize one of the columns in %s for table %s." % \
                (cols, tbl_name)
        if self.sqltype is sqltypes.POSTGRESQL and CAN_COPY:
            conn = self.engine.raw_connection()
            mngr = CopyManager(conn, tbl_name, cols)
            mngr.copy(data, BytesIO)
            conn.commit()
        else:
            #TODO: use bulk insert mappings?
            print("WARNING: You are not using postresql or do not have pgcopy, so this will likely be very slow.")
            self.insert_many('tbl_name', [dict(zip(cols, ro)) for ro in data])


    def filter_query(self, tbls, *args):
        "Query a table and filter results."
        self.grab_session()
        if isiterable(tbls) and not isinstance(tbls, dict):
            if isinstance(tbls[0], type(self.Base)):
                query_args = tbls
            elif isinstance(tbls[0], str):
                query_args = [self.tables[tbl] for tbl in tbls]
            else:
                raise InputError('Unrecognized table specification type: %s.' %
                                 type(tbls[0]))
        else: 
            if isinstance(tbls, type(self.Base)):
                query_args = [tbls]
            elif isinstance(tbls, str):
                query_args = [self.tables[tbls]]
            else:
                raise InputError('Unrecognized table specification type: %s.' %
                                 type(tbls))
                
        return self.session.query(*query_args).filter(*args)


    def select_one(self, tbls, *args):
        """Select the first value that matches requirements given in kwargs 
        from table indicated by tbl_name.
        
        Note that if your specification yields multiple results, this method 
        will just return the first result without exception.
        """
        return self.filter_query(tbls, *args).first()


    def select(self, tbls, *args):
        """Select any and all entries from table given by tbl_name.
        
        The results will be filtered your keyword arguments, as for example
        adding `source='Springer'` to a selection from text content.
        """
        return self.filter_query(tbls, *args).all()


    def insert_db_stmts(self, stmts, db_name):
        "Insert statement, their database, and any affiliated agents."
        # Insert the db info
        print("Adding db %s." % db_name)
        db_ref_id = self.insert(
            'db_info', 
            db_name = db_name, 
            timestamp=get_timestamp()
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
                    assert False, "Unhandled agent role."

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
        "Get abstracts using teh pmids in pmid_list."
        abst_list = self.filter_query(
            [self.TextRef, self.TextContent],
            self.TextContent.text_ref_id==self.TextRef.id,
            self.TextContent.text_type=='abstract',
            self.TextRef.pmid.in_(pmid_list)
            ).all()
        if unzip:
            def unzip_func(s):
                return unzip_string(s.tobytes())
        else:
            def unzip_func(s):
                return s
        return [(r.pmid, unzip_func(c.content)) for (r,c) in abst_list]


    def get_auth_xml_pmcids(self):
        tref_list = self.filter_query(
            [self.TextRef, self.TextContent],
            self.TextRef.id==self.TextContent.text_ref_id,
            self.TextContent.text_type=='fulltext',
            self.TextContent.source=='pmc_auth'
            )
        return [tref.pmcid for tref in tref_list]


    def get_all_pmids(self):
        "Get a list of all the pmids on record."
        return self.get_values(self.select('text_ref'), 'pmid')


def get_aws_db():
    try:
        db = DatabaseManager(DEFAULT_AWS_HOST)
        db.grab_session()
    except Exception as e:
        print("WARNING: Could not connect to aws database due to error:\n%s." % e)
        class DummyManager(object):
            def __getattribute__(self, attr_name, *args, **kwargs):
                print(
                    "WARNING: Unable to get %s because db could not be loaded." % attr_name
                    )
        db = DummyManager()
    return db

'''
def get_abstract_pmids():
    conn = get_connection()
    cur = conn.cursor()
    sql = """SELECT text_ref.pmid FROM text_ref, text_content
                WHERE text_ref.id = text_content.text_ref_id
                AND text_content.text_type = 'abstract';"""
    cur.execute(sql)
    conn.commit()
    return [p[0] for p in cur.fetchall()]


def select_stmts_by_gene(hgnc_name, role=None):
    hgnc_id = hgnc_client.get_hgnc_id(hgnc_name)
    conn = get_connection()
    cur = conn.cursor()
    sql = """SELECT statements.json FROM statements, agents
                WHERE agents.db_name = 'HGNC'
                  AND agents.db_id = %s
                  AND agents.stmt_id = statements.id;"""
    cur.execute(sql, (hgnc_id,))
    conn.commit()
    stmts = []
    for result in cur.fetchall():
        stmt_json_str = result[0]
        stmt_json = json.loads(stmt_json_str)
        stmts.append(Statement._from_json(stmt_json))
    return stmts


def select_text_no_reach():
    """Get text_ref records that have not had any content read by REACH."""
    conn = get_connection()
    cur = conn.cursor()
    sql = """SELECT id, content_type FROM text_content
             WHERE NOT EXISTS
                (SELECT * FROM reach
                 WHERE reach.text_content_id = text_content.id);"""
    cur.execute(sql)
    conn.commit()
    return cur.fetchall()
'''
