import json
import psycopg2 as pg
from indra.statements import *
from indra.databases import hgnc_client

conn = None

def get_connection():
    global conn
    if conn is None:
        aws_host = 'indradb.cwcetxbvbgrf.us-east-1.rds.amazonaws.com'
        conn = pg.connect(host=aws_host, database='indra_db',
                          user='indra_db_user', password='indra_db_pass')
    return conn

def create_tables():
    """Create the tables for the INDRA database."""
    conn = get_connection()
    sql = """
    CREATE TABLE text_ref (
        id serial PRIMARY KEY,
        source VARCHAR NOT NULL,
        pmid VARCHAR UNIQUE,
        pmcid VARCHAR UNIQUE,
        doi VARCHAR UNIQUE,
        url VARCHAR UNIQUE,
        manuscript_id VARCHAR UNIQUE,
        journal VARCHAR,
        publisher VARCHAR,
        pub_date DATE
    );
    CREATE TABLE text_content (
        id serial PRIMARY KEY,
        text_ref_id int4 NOT NULL REFERENCES text_ref(id),
        content_type VARCHAR NOT NULL,
        content TEXT NOT NULL
    );
    CREATE TABLE  db_info (
        id serial PRIMARY KEY,
        db_name VARCHAR NOT NULL,
        timestamp TIMESTAMP DEFAULT now()
    );
    CREATE TABLE statements (
        id serial PRIMARY KEY,
        uuid VARCHAR UNIQUE NOT NULL,
        db_ref int4 REFERENCES db_info(id),
        type VARCHAR NOT NULL,
        json TEXT NOT NULL
    );
    CREATE TABLE agents (
        id serial PRIMARY KEY,
        stmt_id int4 REFERENCES statements(id),
        db_name VARCHAR NOT NULL,
        db_id VARCHAR NOT NULL,
        role VARCHAR NOT NULL
    );
    SET timezone = 'EST'
    """
    cur = conn.cursor()
    cur.execute(sql)
    conn.commit()

def drop_tables():
    """Drop all tables in the INDRA database."""
    conn = get_connection()
    cur = conn.cursor()
    drop_cmds = ['DROP TABLE agents;',
                 'DROP TABLE statements;',
                 'DROP TABLE db_info;',
                 'DROP TABLE text_content;',
                 'DROP TABLE text_ref;',]
    for cmd in drop_cmds:
        try:
            cur.execute(cmd)
            conn.commit()
        except pg.ProgrammingError as pe:
            print(pe)
            conn.rollback()
    conn.commit()

def show_tables():
    """Show all tables in the INDRA database."""
    conn = get_connection()
    cur = conn.cursor()
    cur.execute('SELECT * FROM pg_catalog.pg_tables')
    conn.commit()
    for table_info in cur.fetchall():
        if table_info[2] == 'indra_db_user':
            print(table_info)

def insert_text_ref(**kwargs):
    conn = get_connection()
    # First create the text ref entry
    sql = """INSERT INTO text_ref
                 VALUES (DEFAULT, %(source)s, %(pmid)s, %(pmcid)s, %(doi)s,
                         %(url)s, %(manuscript_id)s, %(journal)s,
                         %(publisher)s, %(pub_date)s)
                 RETURNING id;"""
    cur = conn.cursor()
    args = dict(zip(['source', 'pmid', 'pmcid', 'doi', 'url', 'manuscript_id',
                      'journal', 'publisher', 'pub_date'], [None]*9))
    args.update(kwargs)
    cur.execute(sql, args)
    text_ref_id = cur.fetchone()[0]
    conn.commit()
    return text_ref_id

def insert_text_content(text_ref_id, content_type, content):
    conn = get_connection()
    sql = "INSERT INTO text_content VALUES (DEFAULT, %s, %s, %s);"
    cur = conn.cursor()
    cur.execute(sql, (text_ref_id, content_type, content))
    conn.commit()

def select_all(table):
    conn = get_connection()
    cur = conn.cursor()
    sql = "SELECT * FROM %s;" % table
    cur.execute(sql)
    conn.commit()
    for text_ref in cur.fetchall():
        print(text_ref)

def get_auth_xml_pmcids():
    conn = get_connection()
    cur = conn.cursor()
    sql = """SELECT text_ref.pmcid FROM text_ref, text_content
                WHERE text_ref.id = text_content.text_ref_id
                AND text_content.content_type = 'pmc_auth_xml';"""
    cur.execute(sql)
    conn.commit()
    for pmcid in cur.fetchall():
        print(pmcid)

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
        print("Inserting stmt %s (%d of %d)" % (stmt, stmt_ix+1, len(stmts)))
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



if __name__ == '__main__':
    import pickle
    with open('bel_mapk.pkl', 'rb') as f:
        stmts = pickle.load(f)

