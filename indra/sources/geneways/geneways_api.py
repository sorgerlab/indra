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
import os
from indra.sources.geneways.processor import GenewaysProcessor


# Path to the INDRA data folder
path_this = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(path_this, '../../../data')


def process_geneways_files(input_folder=data_folder, get_evidence=True):
    """Reads in Geneways data and returns a list of statements.

    Parameters
    ----------
    input_folder : Optional[str]
        A folder in which to search for Geneways data. Looks for these
        Geneways extraction data files: human_action.txt,
        human_actionmention.txt, human_symbols.txt.
        Omit this parameter to use the default input folder which is
        indra/data.
    get_evidence : Optional[bool]
        Attempt to find the evidence text for an extraction by downloading
        the corresponding text content and searching for the given offset
        in the text to get the evidence sentence. Default: True

    Returns
    -------
    gp : GenewaysProcessor
        A GenewaysProcessor object which contains a list of INDRA statements
        generated from the Geneways action mentions.
    """
    gp = GenewaysProcessor(input_folder, get_evidence)
    return gp
