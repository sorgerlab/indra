import pickle
import pandas
from indra.sources.evex.processor import get_sentence_for_offset, \
    EvexProcessor, get_sentence_relative_offset
from indra.statements.validate import assert_valid_statements
from . import get_resource_file


def test_process_relations():
    standoff_tar_gz = get_resource_file('evex_test_annots.tar.gz')
    relations_df = pandas.read_csv(get_resource_file('evex_rels.tsv'), sep='\t')
    articles_df = pandas.read_csv(get_resource_file('evex_articles.tsv'),
                                  sep='\t')
    standoff_index = {}
    for aid in articles_df.article_id:
        paper_prefix, paper_id = aid.split(': ')
        key = (
            'pubmed' if paper_prefix == 'PMID' else 'pmc',
            paper_id if paper_prefix == 'PMID' else paper_id.replace('PMC', '')
        )
        standoff_index[key] = standoff_tar_gz

    ep = EvexProcessor(relations_df, articles_df, standoff_index)
    ep.process_statements()
    assert_valid_statements(ep.statements)
    assert len(ep.statements) == 12
    for stmt in ep.statements:
        assert len(stmt.evidence) == 1
        ev = stmt.evidence[0]
        assert ev.text
        assert ev.text_refs
        for agent in stmt.agent_list():
            assert 'EGID' in agent.db_refs
            assert 'TEXT' in agent.db_refs
    return ep


def test_get_sentence_offset():
    text_lines = [('Interferon-gamma regulates alpha 2-macroglobulin and '
                   'alpha 1-antichymotrypsin expression on the '
                   'pretranslational level in HepG2 cells.')]
    line_offsets = [0]

    assert get_sentence_for_offset(text_lines, line_offsets, 17) == \
        text_lines[0]

    text_lines = ['a', 'b', 'c', 'd', 'e', 'f']
    line_offsets = [0, 188, 376, 627, 823, 1129]
    assert get_sentence_for_offset(text_lines, line_offsets, 535) == 'c'


def test_relative_offset():
    line_offsets = [0, 10, 20]
    assert get_sentence_relative_offset(line_offsets, 5) == 5
    assert get_sentence_relative_offset(line_offsets, 15) == 5
    assert get_sentence_relative_offset(line_offsets, 25) == 5


def test_binding_standoff():
    with open(get_resource_file('evex_binding_standoff.pkl'), 'rb') as fh:
        standoff = pickle.load(fh)

    source_id = '19'
    target_id = '20'
    regs = standoff.find_potential_regulations(source_id, target_id)
    assert regs
