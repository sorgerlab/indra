import os
import pickle
import random

from nose.plugins.attrib import attr

from indra.literature import pubmed_client as pubc

from indra.db import util as dbu
from indra.db import client as dbc

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class _PrePaDatabaseTestSetup(object):
    """This object is used to setup the test database into various configs."""
    def __init__(self, max_total_stmts):
        self.test_db = dbu.get_test_db()
        self.test_db._clear(force=True)
        with open(os.path.join(THIS_DIR, 'db_pa_test_input_1M.pkl'), 'rb') as f:
            self.test_data = pickle.load(f)

        if max_total_stmts < len(self.test_data['raw_statements']['tuples']):
            self.stmt_tuples = random.sample(
                self.test_data['raw_statements']['tuples'],
                max_total_stmts
                )
        else:
            self.stmt_tuples = self.test_data['raw_statements']['tuples']

        self.used_stmt_tuples = set()
        return

    def get_available_stmt_tuples(self):
        return list(set(self.stmt_tuples) - self.used_stmt_tuples)

    def load_background(self):
        """Load in all the background provenance metadata (e.g. text_ref).

        Note: This must be done before you try to load any statements.
        """
        # Abbreviate this variable to avoid excessive line breaks.
        td = self.test_data
        tables = ['text_ref', 'text_content', 'reading', 'db_info']

        # Handle the case where we aren't using all the statements.
        if len(self.stmt_tuples) < len(td['raw_statements']['tuples']):
            # Look up the indices for easy access.
            rdg_idx = td['raw_statements']['cols'].index('reading_id')
            tc_idx = td['reading']['cols'].index('text_content_id')
            tr_idx = td['text_content']['cols'].index('text_ref_id')

            # Select only the necessary refs
            inputs = {tbl: set() for tbl in tables}

            # Take all the db_info (there aren't many).
            inputs['db_info'] = set(td['db_info']['tuples'])

            # Filter out un-needed reading provenance.
            for stmt_tpl in self.stmt_tuples:
                rid = stmt_tpl[rdg_idx]
                if not rid:
                    continue
                # Select the reading.
                rdg_tpl = td['reading']['dict'][stmt_tpl[rdg_idx]]
                inputs['reading'].add(rdg_tpl)

                # Select the text content.
                tc_tpl = td['text_content']['dict'][rdg_tpl[tc_idx]]
                inputs['text_content'].add(tc_tpl)

                # Select the text ref.
                inputs['text_ref'].add(td['text_ref']['dict'][tc_tpl[tr_idx]])
        else:
            inputs = {tbl: set(td[tbl]['tuples']) for tbl in tables}

        # Insert the necessary content.
        for tbl in tables:
            print("Loading %s..." % tbl)
            self.test_db.copy(tbl, inputs[tbl], self.test_data[tbl]['cols'])
        return

    def insert_the_statements(self, input_tuples):
        print("Loading %d statements..." % len(input_tuples))
        self.test_db.copy('raw_statements', [t[1:] for t in input_tuples],
                          self.test_data['raw_statements']['cols'][1:])
        print("Inserting agents...")
        dbu.insert_agents(self.test_db, 'raw')
        return

    def add_statements(self):
        """Add statements and agents to the database."""
        input_tuples = self.get_available_stmt_tuples()
        self.insert_the_statements(input_tuples)
        self.used_stmt_tuples |= set(input_tuples)
        return


def _get_prepped_db(num_stmts):
    dts = _PrePaDatabaseTestSetup(num_stmts)
    dts.load_background()
    dts.add_statements()
    return dts.test_db


@attr('nonpublic', 'slow')
def test_get_statements():
    num_stmts = 10000
    db = _get_prepped_db(num_stmts)

    # Test getting all statements
    stmts = dbc.get_statements([], preassembled=False, db=db)
    assert len(stmts) == num_stmts, len(stmts)

    stmts = dbc.get_statements([db.RawStatements.reading_id.isnot(None)],
                               preassembled=False, db=db)
    pmids = {s.evidence[0].pmid for s in random.sample(stmts, 200)}
    assert pmids
    assert None not in pmids
    md_list = pubc.get_metadata_for_ids(list(pmids))
    assert len(md_list) == len(pmids), (len(md_list), len(pmids))

    # Test getting some statements
    stmt_uuid = stmts[0].uuid
    stmts = dbc.get_statements([db.RawStatements.uuid != stmt_uuid],
                               preassembled=False, db=db)
    assert len(stmts) == num_stmts-1, len(stmts)

    # Test getting statements without fix refs.
    stmts = dbc.get_statements([db.RawStatements.reading_id.isnot(None),
                                db.RawStatements.reading_id == db.Reading.id,
                                db.Reading.reader == 'SPARSER'],
                               preassembled=False, fix_refs=False, db=db)
    assert 0 < len(stmts) < num_stmts, len(stmts)
    pmids = {s.evidence[0].pmid for s in random.sample(stmts, 200)}
    assert None in pmids, pmids


@attr('nonpublic', 'slow')
def test_get_statements_by_grot():
    """Test get statements by gene-role-type."""
    num_stmts = 10000
    db = _get_prepped_db(num_stmts)

    stmts = dbc.get_statements_by_gene_role_type('MAP2K1', preassembled=False)
    assert stmts

    stmts = dbc.get_statements_by_gene_role_type('MEK', agent_ns='FPLX',
                                                 preassembled=False)
    assert stmts

    stmts = dbc.get_statements_by_gene_role_type('MAP2K1', preassembled=False,
                                                 fix_refs=False)
    assert stmts

    stmts = dbc.get_statements_by_gene_role_type('MAP2K1', preassembled=False,
                                                 essentials_only=True)
    assert stmts

