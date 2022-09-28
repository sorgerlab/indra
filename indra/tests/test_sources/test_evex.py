import pickle
from indra.sources.evex.processor import get_sentence_for_offset
from . import get_resource_file


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


def test_binding_standoff():
    with open(get_resource_file('evex_binding_standoff.pkl'), 'rb') as fh:
        standoff = pickle.load(fh)

    source_id = '19'
    target_id = '20'
    regs = standoff.find_potential_regulations(source_id, target_id)
    assert regs