import psycopg2 as pg

conn = None

def create_tables():
    """Create the tables for the INDRA database."""
    conn = get_connection()
    sql = """
    CREATE TABLE text_ref (
        id serial PRIMARY KEY,
        source VARCHAR(20) NOT NULL,
        pmid VARCHAR(20) UNIQUE,
        pmcid VARCHAR(20) UNIQUE,
        doi VARCHAR(1000) UNIQUE,
        url VARCHAR(1000) UNIQUE,
        manuscript_id VARCHAR(20) UNIQUE,
        journal VARCHAR(1000),
        publisher VARCHAR(1000),
        pub_date DATE
    );

    CREATE TABLE text_content (
        id serial PRIMARY KEY,
        text_ref_id int4 NOT NULL REFERENCES text_ref(id),
        content_type VARCHAR(100) NOT NULL,
        content TEXT NOT NULL
    );
    """
    cur = conn.cursor()
    cur.execute(sql)
    conn.commit()

def drop_tables():
    """Drop all tables in the INDRA database."""
    conn = get_connection()
    cur = conn.cursor()
    cur.execute('DROP TABLE text_content;')
    cur.execute('DROP TABLE text_ref;')
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


def add_text_ref(**kwargs):
    conn = get_connection()
    # First create the text ref entry
    sql = """INSERT INTO text_ref
                 VALUES (DEFAULT, %(source)s, %(pmid)s, %(pmcid)s, %(doi)s,
                         %(url)s, %(manuscript_id)s, %(journal)s,
                         %(publisher)s, %(pub_date)s);"""
    cur = conn.cursor()
    args = dict(zip(['source', 'pmid', 'pmcid', 'doi', 'url', 'manuscript_id',
                      'journal', 'publisher', 'pub_date'], [None]*9))
    args.update(kwargs)
    cur.execute(sql, args)
    conn.commit()
    print(args)

def add_text_content_by_pmid(pmid, content_type, content):
    conn = get_connection()
    sql = """INSERT INTO text_content
            VALUES (DEFAULT,
                    (SELECT id from text_ref WHERE pmid=%s), %s, %s);"""
    cur = conn.cursor()
    cur.execute(sql, (pmid, content_type, content))
    conn.commit()

def select_all(table):
    conn = get_connection()
    cur = conn.cursor()
    sql = "SELECT * FROM %s;" % table
    cur.execute(sql)
    conn.commit()
    for text_ref in cur.fetchall():
        print(text_ref)

def get_connection():
    global conn
    if conn is None:
        aws_host = 'indradb.cwcetxbvbgrf.us-east-1.rds.amazonaws.com'
        conn = pg.connect(host=aws_host, database='indra_db',
                          user='indra_db_user', password='indra_db_pass')
    return conn

