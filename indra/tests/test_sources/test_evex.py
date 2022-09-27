from indra.sources.evex.processor import get_sentence_for_offset


def test_get_sentence_offset():
    text_lines = [('Interferon-gamma regulates alpha 2-macroglobulin and '
                   'alpha 1-antichymotrypsin expression on the '
                   'pretranslational level in HepG2 cells.')]
    line_offsets = [0]

    assert get_sentence_for_offset(text_lines, line_offsets, 17) == \
        text_lines[0]
