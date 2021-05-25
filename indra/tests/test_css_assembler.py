from indra.assemblers.html.assembler import DEFAULT_SOURCE_COLORS
from indra.assemblers.html.css_assembler import SourceBadgeStyleSheet, \
    StmtsViewStyleSheet, BaseTemplateStyleSheet


def test_source_badge_sheet():
    sbss = SourceBadgeStyleSheet()
    css_str = sbss.make_model()

    # Loop all the colors and test if they are in the generated stylesheet str
    for category, data in DEFAULT_SOURCE_COLORS:
        for source, font_color in data['sources'].items():
            style_str = f".source-{source} " \
                        "{\n" \
                        f"    background-color: {font_color};\n" \
                        f"    color: {data['color']};\n" \
                        "}"
            assert style_str in css_str
