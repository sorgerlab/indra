"""
This module provides a simplified API for invoking the Geneways input processor
, which converts extracted information collected with Geneways into INDRA
statements.

See publication:
Rzhetsky, Andrey, Ivan Iossifov, Tomohiro Koike, Michael Krauthammer, Pauline
Kra, Mitzi Morris, Hong Yu et al. "GeneWays: a system for extracting,
analyzing, visualizing, and integrating molecular pathway data."
Journal of biomedical informatics 37, no. 1 (2004): 43-53.
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.sources.geneways.processor import GenewaysProcessor

def process_geneways(search_path=None):
    """Reads in Geneways data and returns a list of statements.

    Parameters
    ----------
    search_path : list
        a list of directories in which to search for Geneways data.
        Looks for these Geneways extraction data files:
        human_action.txt, human_actionmention.txt,
        human_symbols.txt. Omit this parameter to use the default search path.

    Returns
    -------
    statements : list[indra.statements.Statement]
        A list of INDRA statements generated from the Geneways action mentions.
    """
    if search_path is None:
        search_path = ['./data', '../data', '../../data', '~/data', '.']

    processor = GenewaysProcessor(search_path)
    return processor.statements

