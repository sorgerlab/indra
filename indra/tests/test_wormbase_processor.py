import os
from indra.statements import Complex
from indra.sources.wormbase import process_from_files


this_dir = os.path.dirname(__file__)
test_file_gen = os.path.join(this_dir,
                             'wormbase_tests_data/INTERACTION-GEN_WB_test.tsv')
test_file_mol = os.path.join(this_dir,
                             'wormbase_tests_data/INTERACTION-MOL_WB_3_test.tsv')
test_file_map = os.path.join(this_dir,
                             'wormbase_tests_data/wormbase_entrez_mappings.tsv')


def test_processor():
    processor = process_from_files(test_file_gen, test_file_mol,
                                   test_file_map)
    stmts = processor.statements
    assert len(stmts) == 2
    assert all(isinstance(stmt, Complex) for stmt in stmts)
    for stmt in stmts:
        for agent in stmt.agent_list():
            assert set(agent.db_refs.keys()) >= {'WB', 'EGID'}
