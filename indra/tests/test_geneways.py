from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.sources.geneways.geneways_symbols_parser import \
        GenewaysSymbols
from indra.sources.geneways.geneways_action_parser import \
        GenewaysActionParser
from indra.sources.geneways.geneways_actionmention_parser import \
        GenewaysActionMentionParser
from indra.sources.geneways.geneways_api import process_geneways_files

# Path to the Geneways test/dummy data folder
path_this = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(path_this, 'geneways_tests_data')
symbols_file = os.path.join(data_folder, 'human_symbols.txt')
actionmention_file = os.path.join(data_folder, 'human_actionmention.txt')

def test_geneways_symbols_parser():
    symbols = GenewaysSymbols(symbols_file)

    print(symbols.symbol_to_id('Akt') == ['1'])
    assert(symbols.symbol_to_id('Akt') == ['1'])
    assert(symbols.symbol_to_id('c-Src') == ['2'])

    assert(symbols.id_to_symbol('1') == 'Akt')
    assert(symbols.id_to_symbol('2') == 'c-Src')

    assert(len(symbols.symbols_with_multiple_ids()) == 0)

def test_geneways_actionmention_parser():
    parser = GenewaysActionMentionParser(actionmention_file)

    assert(len(parser.action_mentions) == 4)
    mention0 = parser.action_mentions[0]
    mention1 = parser.action_mentions[1]
    mention2 = parser.action_mentions[2]
    mention3 = parser.action_mentions[3]

    #Make sure that the parser reads in the TSV file into the correct fields

    #First action mention
    assert(mention0.hiid == '1')
    assert(mention0.actionmentionid == '1')
    assert(mention0.negative == '0')
    assert(mention0.upstream == 'c-Src')
    assert(mention0.actiontype == 'phosphorylate')
    assert(mention0.downstream == 'Akt')
    assert(mention0.pmid == '19262695')
    assert(mention0.isFullText == '1')
    assert(mention0.sentencenumber == '4')
    assert(mention0.score == '0.56')
    assert(mention0.prec == '0.78')

    #Second action mention
    assert(mention1.hiid == '1')
    assert(mention1.actionmentionid == '2')
    assert(mention1.negative == '0')
    assert(mention1.upstream == 'c-Src')
    assert(mention1.actiontype == 'phosphorylate')
    assert(mention1.downstream == 'Akt')
    assert(mention1.pmid == '2')
    assert(mention1.isFullText == '0')
    assert(mention1.sentencenumber == '7')
    assert(mention1.score == '0.23')
    assert(mention1.prec == '0.34')

    #Third action mention
    assert(mention2.hiid == '2')
    assert(mention2.actionmentionid == '3')
    assert(mention2.negative == '1')
    assert(mention2.upstream == 'A')
    assert(mention2.actiontype == 'bind')
    assert(mention2.downstream == 'B')
    assert(mention2.pmid == '0')
    assert(mention2.isFullText == '0')
    assert(mention2.sentencenumber == '0')
    assert(mention2.score == '0.12')
    assert(mention2.prec == '0.48')

    #Fourth action mention
    assert(mention3.hiid == '3')
    assert(mention3.actionmentionid == '4')
    assert(mention3.negative == '0')
    assert(mention3.upstream == 'C')
    assert(mention3.actiontype == 'bind')
    assert(mention3.downstream == 'D')
    assert(mention3.pmid == '0')
    assert(mention3.isFullText == '0')
    assert(mention3.sentencenumber == '0')
    assert(mention3.score == '0.22')
    assert(mention3.prec == '0.55')

def test_geneways_action_parser():
    parser = GenewaysActionParser(data_folder)

    actions = parser.actions
    assert(len(actions) == 3)

    action0 = actions[0]
    action1 = actions[1]
    action2 = actions[2]

    assert(action0.hiid == '1')
    assert(action0.up == '2')
    assert(action0.dn == '1')
    assert(action0.actiontype == 'phosphorylate')
    assert(action0.action_count == '1')
    assert(action0.actionmention_count == '2')
    assert(action0.plo == 'P')
    assert(action0.max_score == '0.77')
    assert(action0.max_prec == '0.88')
    assert(len(action0.action_mentions) == 2)

    assert(action1.hiid == '2')
    assert(action1.up == '3')
    assert(action1.dn == '4')
    assert(action1.actiontype == 'bind')
    assert(action1.action_count == '1')
    assert(action1.actionmention_count == '1')
    assert(action1.plo == 'P')
    assert(action1.max_score == '0.12')
    assert(action1.max_prec == '0.34')
    assert(len(action1.action_mentions) == 1)

    assert(action2.hiid == '3')
    assert(action2.up == '5')
    assert(action2.dn == '6')
    assert(action2.actiontype == 'bind')
    assert(action2.action_count == '1')
    assert(action2.actionmention_count == '1')
    assert(action2.plo == 'P')
    assert(action2.max_score == '0.16')
    assert(action2.max_prec == '0.17')
    assert(len(action2.action_mentions) == 1)

def test_geneways_processor():
    processor = process_geneways_files(data_folder, get_evidence=False)

    statements = processor.statements
    assert(len(statements) == 3)

    statement0 = statements[0]
    statement1 = statements[1]
    statement2 = statements[2]

    assert(isinstance(statement0, Phosphorylation))
    assert(statement0.enz.db_refs['TEXT'] == 'c-Src')
    assert(statement0.enz.name == 'A2M')
    assert(statement0.sub.db_refs['TEXT'] == 'Akt')
    assert(statement0.sub.name == 'A1BG')
    assert(statement0.evidence[0].pmid == '19262695')
    assert(statement0.evidence[0].epistemics['direct'])

    assert(isinstance(statement1, Phosphorylation))
    assert(statement1.enz.db_refs['TEXT'] == 'c-Src')
    assert(statement1.enz.name == 'A2M')
    assert(statement1.sub.db_refs['TEXT'] == 'Akt')
    assert(statement1.sub.name == 'A1BG')
    assert(statement1.evidence[0].pmid == '2')
    assert(statement1.evidence[0].epistemics['direct'])

    assert(isinstance(statement2, Complex))
    assert(statement2.members[0].db_refs['TEXT'] == 'C')
    assert(statement2.members[1].db_refs['TEXT'] == 'D')

