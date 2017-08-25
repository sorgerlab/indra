import json
import psycopg2 as pg
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
    try:
        ret = datetime.utcnow().timestamp()
    except:
        now = datetime.utcnow()
        ret = time.mktime(now.timetuple())+now.microsecond/1000000.0
    return ret

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
            UniqueConstraint('pmid', 'doi')
            UniqueConstraint('pmcid', 'doi')
        
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
            UniqueConstraint('text_ref_id', 'source', 'format', 'text_type')
        
        class Reach(self.Base):
            __tablename__ = 'reach'
            id = Column(Integer, primary_key=True)
            text_content_id = Column(Integer, 
                                     ForeignKey('text_content.id'), 
                                     nullable=False)
            text_content = relationship(TextContent)
            version = Column(String(20), nullable=False)
            json = Column(Bytea, nullable=False)
            UniqueConstraint('text_content_id', 'version')
        
        class DBInfo(self.Base):
            __tablename__ = 'db_info'
            id = Column(Integer, primary_key=True)
            db_name = Column(String(20), nullable=False)
            #FIXME: I don't think the default here is really correct.
            timestamp = Column(TIMESTAMP, default=get_timestamp())
        
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

    def insert_text_ref(self, **kwargs):
        "Insert a text_ref entry."
        inputs = dict.fromkeys(self.get_columns('text_ref'))
        inputs.update(kwargs)
        new_text_ref = self.tables['text_ref'](**inputs)
        self.session.add(new_text_ref)
        self.session.commit()
        return new_text_ref.id

    def insert_text_content(self, **kwargs):
        "Insert an entry into text_content"
        inputs = dict.fromkeys(self.get_columns('text_content'))
        inputs.update(kwargs)
        new_text_ref = self.tables['text_content'](**inputs)
        self.session.add(new_text_ref)
        self.session.commit()

    def _filter_query(self, tbl_name, **kwargs):
        "Query a table and filter results."
        args = []
        for attr_name, val in kwargs.items():
            args.append(getattr(self.tables[tbl_name], attr_name) == val)
        return self.session.query(self.tables[tbl_name]).filter(*args)

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
    
    def get_abstracts_by_pmids(self, pmid_list, unzip=True):
        "Get abstracts using teh pmids in pmid_list."
    
    def insert_db_stmts(self, stmts, db_name):
        "Insert statements into the db of db_name"
        
    def get_auth_xml_pmcids(self):
        sql = """SELECT text_ref.pmcid FROM text_ref, text_content
                    WHERE text_ref.id = text_content.text_ref_id AND
                          text_content.text_type = 'fulltext' AND
                          text_content.source = 'pmc_auth';"""
        res = self.session#.something...
        return [p[0] for p in res.all()]
        
    
        
db = DatabaseManager(DEFAULT_AWS_HOST)
db.get_session()

'''
def insert_db_stmts(stmts, db_name):
    # First, add the DB info
    conn = get_connection()
    cur = conn.cursor()
    print("Adding db %s" % db_name)
    sql = """INSERT INTO db_info (db_name) VALUES (%s) RETURNING id;"""
    cur.execute(sql, (db_name,))
    db_ref_id = cur.fetchone()[0]
    # Now, insert the statements
    for stmt_ix, stmt in enumerate(stmts):
S        print("Inserting stmt %s (%d of %d)" % (stmt, stmt_ix+1, len(stmts)))
        sql = """INSERT INTO statements (uuid, db_ref, type, json)
                    VALUES (%s, %s, %s, %s) RETURNING id;"""
        cur.execute(sql, (stmt.uuid, db_ref_id, stmt.__class__.__name__,
                          json.dumps(stmt.to_json())))
        stmt_id = cur.fetchone()[0]
        # Now collect the agents and add them
        for ag_ix, ag in enumerate(stmt.agent_list()):
            if ag is None:
                continue
            if isinstance(stmt, Complex) or \
               isinstance(stmt, SelfModification) or \
               isinstance(stmt, ActiveForm):
                role = 'OTHER'
            elif ag_ix == 0:
                role = 'SUBJECT'
            elif ag_ix == 1:
                role = 'OBJECT'
            else:
                assert False, "Unhandled agent role."
            # If no db_refs for the agent, skip the INSERT that follows
            if not ag.db_refs:
                continue
            sql = """INSERT INTO agents (stmt_id, db_name, db_id, role)
                     VALUES """
            args = []
            sql_list = []
            for db, db_id in ag.db_refs.items():
                args += [stmt_id, db, db_id, role]
                sql_list.append('(%s, %s, %s, %s)')
            sql = sql + ','.join(sql_list) + ';'
            cur.execute(sql, args)
    conn.commit()

def select(table, field=None):
    conn = get_connection()
    cur = conn.cursor()
    if field is None:
        # WARNING: SQL injection!
        sql = "SELECT * FROM %s;" % table
        cur.execute(sql)
    else:
        # WARNING: SQL injection!
        sql = "SELECT " + field + " FROM " + str(table) + ";"
        cur.execute(sql)
    conn.commit()
    return cur.fetchall()


def get_text_refs_by_pmcid(pmcid_list):
    conn = get_connection()
    cur = conn.cursor()
    cur.execute("SELECT pmcid, id FROM text_ref WHERE pmcid IN %s;",
                (pmcid_list,))
    conn.commit()
    return cur.fetchall()


def get_abstracts_by_pmids(pmid_list, unzip=True):
    conn = get_connection()
    cur = conn.cursor()
    cur.execute("""SELECT text_ref.pmid, text_content.content
                   FROM text_ref, text_content
                   WHERE text_content.text_type = 'abstract' AND
                         text_content.text_ref_id = text_ref.id AND
                         text_ref.pmid IN %s;""", (tuple(pmid_list),))
    conn.commit()
    results = []
    for pmid, abstract_mem in cur.fetchall():
        results.append((pmid, unzip_string(abstract_mem.tobytes())))
    return results


def get_all_pmids():
    conn = get_connection()
    cur = conn.cursor()
    cur.execute("SELECT pmid FROM text_ref;")
    conn.commit()
    return [r[0] for r in cur.fetchall()]


def get_text_ref_by_id(id_str, id_type):
    """Look up the text ref using an id
    
    Note: This should eventually support looking up by multiple ids.
    """
    conn = get_connection()
    cur = conn.cursor()
    fmt = "SELECT id FROM text_ref WHERE {id_type} = {id}"
    cur.execute(fmt.format(id = id_str, id_type = id_type))
    conn.commit()
    res = cur.fetchone()
    if res:
        return res[0]
    else:
        return None


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




def insert_reach(text_content_id, version, json_str):
    """Insert REACH reading results."""
    conn = get_connection()
    cur = conn.cursor()
    sql = """INSERT INTO reach (text_content_id, version, json)
             VALUES (%s, %s, %s);"""
    cur.execute(sql, (text_content_id, version, json_str))
    conn.commit()


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