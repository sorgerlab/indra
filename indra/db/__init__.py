from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import json
from indra.statements import *
from indra.util import unzip_string
from indra.databases import hgnc_client
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, UniqueConstraint, ForeignKey,\
    TIMESTAMP, create_engine, inspect
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.dialects.postgresql import BYTEA
from datetime import datetime
import time
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


class DatabaseManager(object):
    def __init__(self, host, sqltype='postgresql'):
        self.host = host
        self.Base = declarative_base()
        
        if sqltype is 'postgresql':
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
        self.engine = create_engine(host)
        self.session = None


    def get_connection(self):
        "For backwards compatability."
        self.get_session()


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


    def get_session(self):
        "Get an active session with the database."
        if self.session is None:
            DBSession = sessionmaker(bind=self.engine)
            self.session = DBSession()


    def get_tables(self):
        "Get a list of available tables."
        return [tbl_name for tbl_name in self.tables.keys()]


    def show_tables(self):
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


    def get_values(self, entry_list, col_names):
        "Get the column values from the entries in entry_list"
        ret = []
        for entry in entry_list:
            if isiterable(col_names):
                ret.append([getattr(entry, a) for a in col_names])
            else:
                ret.append(getattr(entry, col_names))
        return ret 


    def insert(self, tbl_name, ret_vals = 'id', **input_dict):
        "Insert a an entry into specified table, and return id."
        inputs = dict.fromkeys(self.get_columns(tbl_name))
        inputs.update(input_dict)
        new_entry = self.tables[tbl_name](**inputs)
        self.session.add(new_entry)
        self.commit("Excepted while trying to insert %s into %s" % 
                  (inputs, tbl_name))
        return self.get_values([new_entry], ret_vals)[0]


    def insert_many(self, tbl_name, input_dict_list, ret_info = 'id'):
        "Insert many records into the table given by table_name."
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


    def insert_text_ref(self, **kwargs):
        "Insert a text_ref entry."
        return self.insert('text_ref', **kwargs)


    def insert_text_content(self, **kwargs):
        "Insert an entry into text_content"
        return self.insert('text_content', **kwargs)


    def insert_db_stmts(self, stmts, db_name):
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
                if isinstance(stmt, Complex) or \
                    isinstance(stmt, SelfModification) or \
                    isinstance(stmt, ActiveForm):
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


    def _filter_query(self, tbl_name, **kwargs):
        """Query a table and filter results.
        
        Note: it is assumed, currently, that iterable arguments are intended
        for use with an `in` clause, which is reasonable given the current
        schema, however, if that should change, this may need to change.
        """
        args = []
        for attr_name, val in kwargs.items():
            attr = getattr(self.tables[tbl_name], attr_name)
            if isiterable(val):
                args.append(attr.in_(val))
            else:
                args.append(attr == val)
        return self.session.query(self.tables[tbl_name]).filter(*args)


    def select_one(self, tbl_name, **kwargs):
        """Select the first value that matches requirements given in kwargs 
        from table indicated by tbl_name.
        
        Note that if your specification yields multiple results, this method 
        will just return the first result without exception.
        """
        return self._filter_query(tbl_name, **kwargs).first()


    def select(self, tbl_name, **kwargs):
        """Select any and all entries from table given by tbl_name.
        
        The results will be filtered your keyword arguments, as for example
        adding `source='Springer'` to a selection from text content.
        """
        return self._filter_query(tbl_name, **kwargs).all()


    def get_text_ref_by_pmcid(self, pmcid):
        "Read a text ref using it's pmcid"
        return self._filter_query('text_ref', pmcid=pmcid).first()


    def get_text_refs_by_pmid(self, pmid):
        "Get all the text refs with a given pmid."
        return self._filter_query('text_ref', pmid=pmid).all()


    def get_text_ref_by_id(self, tr_id):
        "Get a text_ref using the text_ref_id."
        return self._filter_query('text_ref', id=tr_id).first()


    def get_text_refs_by_pmcid(self, pmcid_list):
        "Get multipe refs indicated by pmcids in pmcid_list"
        return self._filter_query('text_ref', pmcid=pmcid_list).all()


    def get_abstracts_by_pmids(self, pmid_list, unzip=True):
        "Get abstracts using teh pmids in pmid_list."
        Ref = self.tables['text_ref']
        Cont = self.tables['text_content']
        q = self.session.query(Ref, Cont)
        abst_list = q.filter(
            Ref.pmid.in_(pmid_list), 
            Cont.text_ref_id==Ref.id, 
            Cont.text_type=='abstract'
            ).all()
        if unzip:
            def unzip_func(s):
                return unzip_string(s.tobytes())
        else:
            def unzip_func(s):
                return s
        return [(r.pmid, unzip_func(c.content)) for (r,c) in abst_list]


    def get_auth_xml_pmcids(self):
        sql = """SELECT text_ref.pmcid FROM text_ref, text_content
                    WHERE text_ref.id = text_content.text_ref_id AND
                          text_content.text_type = 'fulltext' AND
                          text_content.source = 'pmc_auth';"""
        res = self.session#.something...
        return [p[0] for p in res.all()]


    def get_all_pmids(self):
        "Get a list of all the pmids on record."
        return self.get_values(self.select('text_ref'), 'pmid')



try:
    db = DatabaseManager(DEFAULT_AWS_HOST)
    db.get_session()
except Exception as e:
    print("WARNING: Could not connect to aws database due to error:\n%s." % e)
    class DummyManager(object):
        def __getattribute__(self, attr_name, *args, **kwargs):
            print(
                "WARNING: Unable to get %s because db could not be loaded." % attr_name
                )
    db = DummyManager()

'''
def get_auth_xml_pmcids():
    conn = get_connection()
    cur = conn.cursor()
    sql = """SELECT text_ref.pmcid FROM text_ref, text_content
                WHERE text_ref.id = text_content.text_ref_id AND
                      text_content.text_type = 'fulltext' AND
                      text_content.source = 'pmc_auth';"""
    cur.execute(sql)
    conn.commit()
    return [p[0] for p in cur.fetchall()]


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
