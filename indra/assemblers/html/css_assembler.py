import logging
from typing import Optional, List, Tuple, Dict, Union
from os.path import dirname, abspath, join

import jinja2
from jinja2 import Environment, FileSystemLoader

from .assembler import make_source_colors, DEFAULT_SOURCE_COLORS

logger = logging.getLogger(__name__)
HERE = dirname(abspath(__file__))

loader = FileSystemLoader(join(HERE, 'templates'))
env = Environment(loader=loader)

sources_template_file = 'indra/source_style_template.css'

SourceColors = List[Tuple[str, Dict[str, Union[str, Dict[str, str]]]]]


class CSSAssembler:
    """Create a stylesheet from a jinja2 template"""
    template_file = NotImplemented

    def __init__(self):
        self.template: jinja2.Template = env.get_template(self.template_file)
        self.model: Optional[str] = None  # Stylesheet as string
        self.template_kwargs: Dict = {}

    def make_model(self) -> str:
        """Render the template"""
        if self.model is None:
            self.model = self.template.render(**self.template_kwargs)
        return self.model

    def save_model(self, fname: str):
        """Save the CSS stylesheet to a file"""
        model = self.make_model()

        with open(fname, 'w') as fh:
            fh.write(model)


class SourceBadgeStyles(CSSAssembler):
    """Stylesheet defining color, background-color for source count badges"""
    template_file = sources_template_file

    def __init__(self, source_colors: SourceColors = DEFAULT_SOURCE_COLORS):
        super().__init__()
        self.template_kwargs = {'source_colors': source_colors}


class StyleSheet:
    """Assemble a stylesheet from a list of CSSAssembler instances"""
    def __init__(self, sheets: List[CSSAssembler]):
        self.sheets: List[CSSAssembler] = sheets
        self.stylesheet: str = ''

    def make_sheet(self) -> str:
        """Join together the output of the provided CSS template assemblers

        Returns
        -------
        str
        """
        if len(self.stylesheet) < 1:
            self.stylesheet = '\n'.join(sheet.make_model() for sheet in
                                        self.sheets)

        return self.stylesheet

    def save_sheet(self, fname: str):
        """Assembles and then saves the stylesheet to the provided file path

        Parameters
        ----------
        fname :
            The filepath to the CSS output file

        """
        sheet = self.make_sheet()

        with open(fname, 'w') as fh:
            fh.write(sheet)
