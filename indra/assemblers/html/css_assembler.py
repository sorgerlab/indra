"""
Assemble CSS rules to be inserted in an HTML file's style tag or to be saved
in a separate CSS file.
"""

import logging
from typing import Optional, Dict
from os.path import dirname, abspath, join

from jinja2 import Environment, FileSystemLoader, Template

from .assembler import DEFAULT_SOURCE_COLORS

logger = logging.getLogger(__name__)
HERE = dirname(abspath(__file__))


# Note that if templates of the same exist in both 'templates/indra' and
# 'templates/macros', the template in the first directory will be used.
loader = FileSystemLoader([join(HERE, 'templates', 'indra'),
                           join(HERE, 'templates', 'macros')])
env = Environment(loader=loader)

stylesheet_template_file = 'standalone_stylesheet_template.css'


class CSSAssembler:
    """Create a stylesheet from a jinja2 template"""
    template_file = stylesheet_template_file

    def __init__(self):
        self.template: Template = env.get_template(self.template_file)
        self.model: Optional[str] = None  # Stylesheet as string
        self.template_kwargs: Dict = NotImplemented

    def make_model(self) -> str:
        """Render the template"""
        if self.model is None:
            self.model = self.template.render(**self.template_kwargs).strip()
        return self.model

    def save_model(self, fname: str):
        """Save the CSS rules to a file"""
        model = self.make_model()

        with open(fname, 'w') as fh:
            fh.write(model)


class SourceBadgeStyleSheet(CSSAssembler):
    """Stylesheet defining color, background-color for source count badges"""

    def __init__(self, source_colors=None):
        super().__init__()
        if source_colors is None:
            source_colors = DEFAULT_SOURCE_COLORS
        self.template_kwargs = {'source_colors': source_colors,
                                'simple': False,
                                'only_source_badges': True,
                                'base_template': False,
                                'stmts_view': False}


class BaseTemplateStyleSheet(CSSAssembler):
    """Stylesheet for the base template"""

    def __init__(self, simple: bool, source_colors=None):
        super().__init__()
        if source_colors is None:
            source_colors = DEFAULT_SOURCE_COLORS
        self.template_kwargs = {'source_colors': source_colors,
                                'simple': simple,
                                'only_source_badges': False,
                                'base_template': True,
                                'stmts_view': False}


class StmtsViewStyleSheet(CSSAssembler):
    """Stylesheet for the statements view template"""

    def __init__(self, simple: bool, source_colors=None):
        super().__init__()
        if source_colors is None:
            source_colors = DEFAULT_SOURCE_COLORS
        self.template_kwargs = {'source_colors': source_colors,
                                'simple': simple,
                                'only_source_badges': False,
                                'base_template': False,
                                'stmts_view': True}
