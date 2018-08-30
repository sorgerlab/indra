import os
import pickle
import random

from nose.plugins.attrib import attr

from indra.literature import pubmed_client as pubc

from indra.db import util as dbu
from indra.db import client as dbc
from indra.statements import stmts_from_json

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

    stmts = dbc.get_statements_by_gene_role_type('MAP2K1', preassembled=False,
                                                 db=db)
    assert stmts

    stmts = dbc.get_statements_by_gene_role_type('MEK', agent_ns='FPLX',
                                                 preassembled=False, db=db)
    assert stmts

    stmts = dbc.get_statements_by_gene_role_type('MAP2K1', preassembled=False,
                                                 fix_refs=False, db=db)
    assert stmts

    stmts = dbc.get_statements_by_gene_role_type('MAP2K1', preassembled=False,
                                                 essentials_only=True, db=db)
    assert stmts


@attr('nonpublic')
def test_get_statement_jsons_by_agent():
    # Note that this deliberately uses the primary (production) database in
    # testing. This is only allowed because only retrieval is tested, however
    # PLEASE PROCEED WITH CARE WHEN MODIFYING THIS TEST.

    # TODO: don't rely on the primary database, because that's scary in a test.
    agents = [(None, 'MEK', 'FPLX'), (None, 'ERK', 'FPLX')]
    stmt_jsons = dbc.get_statement_jsons_from_agents(agents=agents,
                                                     stmt_type='Phosphorylation')
    assert stmt_jsons
    assert stmt_jsons['statements']
    assert stmt_jsons['total_evidence']
    assert stmt_jsons['evidence_returned']
    stmts = stmts_from_json(stmt_jsons['statements'].values())
    assert len(stmts) == len(stmt_jsons['statements'])
    for s in stmts:
        s_agents = [(None, ag_id, ag_ns) for ag in s.agent_list()
                    for ag_ns, ag_id in ag.db_refs.items()]
        for ag_tpl in agents:
            assert ag_tpl in s_agents


def test_get_statement_jsons_options():
    # Test all possible options regarding the number of statements returned.
    # Note that this suffices to test the same options in other related
    # functions as well (e.g. the paper version).
    options = {'max_stmts': 2, 'ev_limit': 4, 'offset': 5, 'best_first': False}
    agents = [('SUBJECT', 'MEK', 'FPLX'), ('OBJECT', 'ERK', 'FPLX')]
    option_dicts = [{}]
    for key, value in options.items():
        nd = {key: value}
        new_option_dicts = []
        for option_dict in option_dicts:
            new_option_dicts.append(option_dict)
            new_option_dicts.append({**option_dict, **nd})
        option_dicts = new_option_dicts

    evidence_count_record = {}
    for option_dict in option_dicts:
        res = dbc.get_statement_jsons_from_agents(agents=agents,
                                                  stmt_type='Phosphorylation',
                                                  **option_dict)
        assert res
        assert len(res['statements'])
        stmts = res['statements']
        if 'max_stmts' in option_dict.keys():
            assert len(stmts) == option_dict['max_stmts']
        else:
            assert len(stmts) == 3

        if 'ev_limit' in option_dict.keys():
            assert all([len(s.evidence) <= options['ev_limit'] for s in stmts.values()])
        else:
            if evidence_count_record:
                assert all([len(s.evidence) == evidence_count_record[mk_hash]
                            for mk_hash, s in stmts.items()])
            else:
                for mk_hash, stmt in stmts.items():
                    evidence_count_record[mk_hash] = len(stmt.evidence)
    return


@attr('nonpublic')
def test_get_statement_jsons_by_paper_id():
    paper_refs = [
        [('pmid', '27769048'), ('pmcid', 'PMC5363599')],
        [('doi', '10.3389/FIMMU.2017.00781')],
        [('pmcid', 'PMC4789553')]
        ]
    stmt_jsons = dbc.get_statement_jsons_from_papers(paper_refs)
    assert stmt_jsons
    assert stmt_jsons['statements']
    assert stmt_jsons['total_evidence']
    stmts = stmts_from_json(stmt_jsons['statements'].values())
    assert len(stmts) == len(stmt_jsons['statements'])
    pmid_set = {ev.pmid for s in stmts for ev in s.evidence}
    assert len(pmid_set) >= len(paper_refs)


def test_get_statement_jsons_by_mk_hash():
    mk_hashes = {-35990550780621697, -34509352007749723, -33762223064440060,
                 -33265410753427801, -33264422871226821, -33006503639209361,
                 -32655830663272427, -32156860839881910, -31266440463655983,
                 -30976459682454095, -30134498128794870, -28378918778758360,
                 -24358695784465547, -24150179679440010, -23629903237028340,
                 -23464686784015252, -23180557592374280, -22224931284906368,
                 -21436209384793545, -20730671023219399, -20628745469039833,
                 -19678219086449831, -19263047787836948, -19233240978956273,
                 -18854520239423344, -18777221295617488, -18371306000768702,
                 -17790680150174704, -17652929178146873, -17157963869106438,
                 -17130129999418301, -16284802852951883, -16037105293100219,
                 -15490761426613291, -14975140226841255, -14082507581435438,
                 -13857723711006775, -12377086298836870, -11313819223154032,
                 -11213416806465629, -10533303510589718, -9966418144787259,
                 -9862339997617041,  -9169838304025767, -7914540609802583,
                 -5761487437008515, -5484899507470794, -4221831802677562,
                 -3843980816183311, -3444432161721189, -2550187846777281,
                 -1690192884583623, -1574988790414009, -776020752709166,
                 -693617835322587, -616115799439746, 58075179102507,
                 1218693303789519, 1858833736757788, 1865941926926838,
                 1891718725870829, 3185457948420843, 3600108659508886,
                 3858621152710053, 4594557398132265, 5499056407723241,
                 6796567607348165, 6828272448940477, 6929632245307987,
                 7584487035784255, 8424911311360927, 8837984832930769,
                 10511090751198119, 10789407105295331, 10924988153490362,
                 11707113199128693, 12528041861567565, 13094138872844955,
                 13166641722496149, 13330125910711684, 13347703252882432,
                 15002261599485956, 16397210433817325, 16975780060710533,
                 17332680694583377, 17888579535249950, 19337587406307012,
                 22774500444258387, 23665225082661845, 23783937267011041,
                 24050979216195140, 24765024299377586, 25290573037450021,
                 29491428193112311, 30289509021065753, 30992174235867673,
                 31766667918079590, 31904387104764159, 34782800852366343,
                 35686927318045812}
    stmt_jsons = dbc.get_statement_jsons_from_hashes(mk_hashes)
    assert stmt_jsons
    assert stmt_jsons['statements']
    assert stmt_jsons['total_evidence']
    stmts = stmts_from_json(stmt_jsons['statements'].values())
    assert len(stmts) == len(stmt_jsons['statements'])
    assert len(stmts) == len(mk_hashes)
