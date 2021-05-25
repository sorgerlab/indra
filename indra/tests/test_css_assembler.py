from typing import List, Tuple, Dict, Any

from indra.assemblers.html.assembler import DEFAULT_SOURCE_COLORS
from indra.assemblers.html.css_assembler import SourceBadgeStyleSheet, \
    StmtsViewStyleSheet, BaseTemplateStyleSheet


def _source_badge_in_str(css_str: str,
                         source_colors: List[Tuple[str, Dict[str, Any]]]) ->\
        bool:
    # Loop all the colors and test if they are in the generated stylesheet str
    for category, data in source_colors:
        for source, font_color in data['sources'].items():
            style_str = f".source-{source} " \
                        "{\n" \
                        f"    background-color: {font_color};\n" \
                        f"    color: {data['color']};\n" \
                        "}"
            assert style_str in css_str
    return True


def test_source_badge_sheet():
    sbss = SourceBadgeStyleSheet()
    css_str = sbss.make_model()

    assert _source_badge_in_str(css_str, DEFAULT_SOURCE_COLORS)


def test_base_template_sheet():
    btss_simple = BaseTemplateStyleSheet(simple=False)
    css_str = btss_simple.make_model()
    assert _source_badge_in_str(css_str, DEFAULT_SOURCE_COLORS)
