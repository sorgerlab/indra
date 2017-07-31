import psycopg2 as pg

def create_tables(conn):
    """Create the tables for the INDRA database."""
    sql = """
    CREATE TABLE text_ref (
        id serial PRIMARY KEY,
        source VARCHAR(20) not null,
        pmid VARCHAR(20),
        pmcid VARCHAR(20),
        doi VARCHAR(1000),
        url VARCHAR(1000),
        manuscript_id VARCHAR(20),
        journal VARCHAR(1000),
        publisher VARCHAR(1000),
        pub_date DATE
    );

    CREATE TABLE text_content (
        id serial PRIMARY KEY,
        text_ref_id int4 REFERENCES text_ref(id),
        content_type VARCHAR(100),
        content TEXT
    );
    """
    cur = conn.cursor()
    cur.execute(sql)
    conn.commit()

def drop_tables(conn):
    """Drop all tables in the INDRA database."""
    cur = conn.cursor()
    cur.execute('DROP TABLE text_content;')
    cur.execute('DROP TABLE text_ref;')
    conn.commit()

def show_tables(conn):
    """Show all tables in the INDRA database."""
    cur = conn.cursor()
    cur.execute('SELECT * FROM pg_catalog.pg_tables')
    conn.commit()
    for table_info in cur.fetchall():
        if table_info[2] == 'indra_db_user':
            print(table_info)

if __name__ == '__main__':
    # Get the connection
    aws_host = 'indradb.cwcetxbvbgrf.us-east-1.rds.amazonaws.com'
    conn = pg.connect(host=aws_host, database='indra_db', user='indra_db_user',
                            password='indra_db_pass')

    # Create the database tables
    #create_tables(conn)
    show_tables(conn)

