from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.sources.tees.tees_api import tees_process_text
from indra.sources.tees.parse_tees import tees_parse_networkx_to_dot
from indra.statements import Phosphorylation, Dephosphorylation

def test_process_phosphorylation():
    """Test the extraction of phosphorylation with a simple example."""
    s = 'Ras leads to the phosphorylation of Braf.'
    tp = tees_process_text(s)
    statements = tp.statements

    #print('Number of extracted statements:', len(statements))

    assert(len(statements) == 1)

    statement = statements[0]
    #print(statement)
    assert(isinstance(statement, Phosphorylation))

    enz = statement.enz
    #print(enz)
    assert(enz.db_refs['TEXT'] == 'Ras')

    sub = statement.sub
    #print(sub)
    assert(sub.db_refs['TEXT'] == 'Braf')

def test_process_dephosphorylation():
    """Test the extraction of a dephosphorylation sentence. This sentence is
    processed into two INDRA statements, at least one of which is correct."""
    s = 'Here we show that febuxostat suppresses LPS-induced MCP-1 production and mRNA expression via activating MAPK phosphatase-1 (MKP-1) which, in turn, leads to dephosphorylation and inactivation of JNK in macrophages.'
    tp = tees_process_text(s)
    statements = tp.statements
    #print('Statements: ', statements)
    tees_parse_networkx_to_dot(tp.G, 'moo.dot', tp.G.node)

    #print('Number of extracted statements:', len(statements))

    # We'll set this to true if at least one of the extracted statements meets
    # our criteria for correctness
    some_statement_correct = False

    # Extracting this particular sentence with TEES produces a couple
    # statements, at least one of which is correct
    for statement in statements:
        statement_correct = False
        enz = statement.enz.db_refs['TEXT']
        sub = statement.sub.db_refs['TEXT']
        #print('ENZYME:', enz)
        #print('SUBSTRATE:', sub)

        entities_correct = enz == 'MAPK phosphatase-1' and sub == 'JNK'
        type_correct = isinstance(statement, Dephosphorylation)
        statement_correct = type_correct and entities_correct
        some_statement_correct = statement_correct or some_statement_correct

    assert(some_statement_correct)

