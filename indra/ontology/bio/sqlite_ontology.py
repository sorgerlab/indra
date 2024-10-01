"""This module implements an SQLite back end to the
INDRA BioOntology."""
import os
import sqlite3
import logging
from collections import defaultdict
from indra.ontology.ontology_graph import IndraOntology
from indra.ontology.bio.ontology import CACHE_DIR
from indra.ontology.bio import bio_ontology


logger = logging.getLogger(__name__)


DEFAULT_SQLITE_ONTOLOGY = os.path.join(CACHE_DIR, 'bio_ontology.db')


class SqliteOntology(IndraOntology):
    def __init__(self, db_path=DEFAULT_SQLITE_ONTOLOGY):
        super().__init__()
        self.db_path = db_path
        build_sqlite_ontology(db_path)
        conn = sqlite3.connect(db_path)
        self.cur = conn.cursor()

    def initialize(self):
        self._initialized = True

    def isa_or_partof(self, ns1, id1, ns2, id2):
        q = """SELECT 1 FROM relationships
               WHERE child_id=? AND child_ns=? AND parent_id=? AND parent_ns=?
               LIMIT 1;"""
        self.cur.execute(q, (id1, ns1, id2, ns2))
        return self.cur.fetchone() is not None

    def child_rel(self, ns, id, rel_types):
        q = """SELECT children FROM child_lookup
               WHERE parent_id=? AND parent_ns=?
               LIMIT 1;"""
        self.cur.execute(q, (id, ns))
        res = self.cur.fetchone()
        if res is None:
            yield from []
        else:
            yield from [tuple(x.split(':', 1)) for x in res[0].split(',')]

    def get_parents(self, ns, id):
        return list(self.parent_rel(ns, id, {'isa', 'partof'}))

    def get_children(self, ns, id):
        return list(self.child_rel(ns, id, {'isa', 'partof'}))

    def parent_rel(self, ns, id, rel_types):
        q = """SELECT parents FROM parent_lookup
               WHERE child_id=? AND child_ns=?
               LIMIT 1;"""
        self.cur.execute(q, (id, ns))
        res = self.cur.fetchone()
        if res is None:
            yield from []
        else:
            yield from [tuple(x.split(':', 1)) for x in res[0].split(',')]

    def get_node_property(self, ns, id, property):
        return None

    def get_id_from_name(self, ns, name):
        return None


def build_sqlite_ontology(db_path=DEFAULT_SQLITE_ONTOLOGY, force=False):
    # If the database already exists and we are not forcing a rebuild, return
    if os.path.exists(db_path) and not force:
        return

    # Initialize the bio ontology and build the transitive closure
    bio_ontology.initialize()
    bio_ontology._build_transitive_closure()

    # Set up connection
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    logger.info('Building SQLite ontology at %s' % db_path)
    # First, we create the relationships table and populate
    # it with child/parent pairs
    q = """CREATE TABLE relationships (
        child_id TEXT NOT NULL,
        child_ns TEXT NOT NULL,
        parent_id TEXT NOT NULL,
        parent_ns TEXT NOT NULL,
        UNIQUE (child_id, child_ns, parent_id, parent_ns)
    );"""
    cur.execute(q)

    # Insert into the database in chunks
    chunk_size = 10000
    tc = sorted(bio_ontology.transitive_closure)
    all_children = defaultdict(set)
    all_parents = defaultdict(set)
    for i in range(0, len(tc), chunk_size):
        chunk = tc[i:i+chunk_size]
        chunk_values = [(child.split(':', 1)[1], child.split(':')[0],
                         parent.split(':', 1)[1], parent.split(':')[0])
                        for child, parent in chunk]
        for cid, cns, pid, pns in chunk_values:
            all_children[(pid, pns)].add('%s:%s' % (cns, cid))
            all_parents[(cid, cns)].add('%s:%s' % (pns, pid))
        cur.executemany("""INSERT INTO relationships (child_id, 
                        child_ns, parent_id, parent_ns) 
                        VALUES (?, ?, ?, ?);""", chunk_values)
    q = """CREATE INDEX idx_child_parent ON relationships 
        (child_id, child_ns, parent_id, parent_ns);"""
    cur.execute(q)

    # Next, create child and parent lookup tables and populate them
    q = """CREATE TABLE child_lookup (
        parent_id TEXT NOT NULL,
        parent_ns TEXT NOT NULL,
        children TEXT NOT NULL,
        UNIQUE (parent_id, parent_ns)
    );"""
    cur.execute(q)
    q = """CREATE TABLE parent_lookup (
        child_id TEXT NOT NULL,
        child_ns TEXT NOT NULL,
        parents TEXT NOT NULL,
        UNIQUE (child_id, child_ns)
    );"""
    cur.execute(q)
    for (pid, pns), children in all_children.items():
        cur.execute("INSERT INTO child_lookup (parent_id, parent_ns, children) "
                    "VALUES (?, ?, ?);",
                    (pid, pns, ','.join(children)))
    for (cid, cns), parents in all_parents.items():
        cur.execute("INSERT INTO parent_lookup (child_id, child_ns, parents) "
                    "VALUES (?, ?, ?);",
                    (cid, cns, ','.join(parents)))
    # Now add indices to the lookup tables
    q = """CREATE INDEX idx_child_lookup ON child_lookup 
        (parent_id, parent_ns);"""
    cur.execute(q)
    q = """CREATE INDEX idx_parent_lookup ON parent_lookup 
        (child_id, child_ns);"""
    cur.execute(q)

    conn.commit()
    conn.close()
    logger.info('Finished building SQLite ontology')
