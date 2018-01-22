from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.sources.geneways.geneways_symbols_parser import \
        GenewaysSymbols
#from indra.sources.geneways.geneways_action_parser import \
#        GenewaysActionParser
#from indra.sources.geneways.geneways_action_mention_parser import \
#        GenewaysActionMentionParser
#from indra.sources.geneways.geneways_api import process_geneways_files
#from indra.sources.geneways.processor import GenewaysProcessor

# Path to the Geneways test/dummy data folder
path_this = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(path_this,
        '../../data/tests_data/geneways_tests_data')
symbols_file = os.path.join(data_folder, 'human_symbols.txt')
actionmention_file = os.path.join(data_folder,
        'human_actionmention.txt')
action_file = os.path.join(data_folder,
        'human_action.txt')

def test_geneways_symbols_parser():
    symbols = GenewaysSymbols(symbols_file)

    print(symbols.symbol_to_id('Akt') == ['1'])
    assert(symbols.symbol_to_id('Akt') == ['1'])
    assert(symbols.symbol_to_id('c-Src') == ['2'])

    assert(symbols.id_to_symbol('1') == 'Akt')
    assert(symbols.id_to_symbol('2') == 'c-Src')

    assert(len(symbols.symbols_with_multiple_ids()) == 0)
