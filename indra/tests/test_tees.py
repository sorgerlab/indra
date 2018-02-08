from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.sources.tees.tees_api import tees_process_text
from indra.sources.tees.parse_tees import tees_parse_networkx_to_dot
from indra.statements import Phosphorylation, Dephosphorylation, \
        IncreaseAmount, DecreaseAmount, Complex

def test_process_phosphorylation():
    #Test the extraction of phosphorylation with a simple example.
    s = 'Ras leads to the phosphorylation of Braf.'
    tp = tees_process_text(s)
    statements = tp.statements
    
    # Should only extract one statement
    assert(len(statements) == 1)

    # The statement should be a Phosphorylation statement
    statement = statements[0]
    assert(isinstance(statement, Phosphorylation))

    # Check that the enzyme and substrate of the statement match the text
    enz = statement.enz
    assert(enz.db_refs['TEXT'] == 'Ras')
    #
    sub = statement.sub
    assert(sub.db_refs['TEXT'] == 'Braf')

def test_process_dephosphorylation():
    #Test the extraction of a dephosphorylation sentence. This sentence is
    #processed into two INDRA statements, at least one of which is correct.

    # Text statement describing phosphorylation
    s = 'Here we show that febuxostat suppresses LPS-induced MCP-1 ' + \
            'production and mRNA expression via activating ' + \
            'MAPK phosphatase-1 (MKP-1) which, in turn, leads to ' + \
            'dephosphorylation and inactivation of JNK in macrophages.'
    tp = tees_process_text(s)
    statements = tp.statements

    # We'll set this to true if at least one of the extracted statements meets
    # our criteria for correctness
    some_statement_correct = False

    # Extracting this particular sentence with TEES produces a couple
    # statements, at least one of which is correct
    for statement in statements:
        statement_correct = False
        enz = statement.enz.db_refs['TEXT']
        sub = statement.sub.db_refs['TEXT']

        # Does this statement have the entities named in the text?
        entities_correct = enz == 'MAPK phosphatase-1' and sub == 'JNK'
        
        # Is the statement a phosphorylation statement?
        type_correct = isinstance(statement, Dephosphorylation)

        # If yes to both, the correct statement was among those extracted
        statement_correct = type_correct and entities_correct
        some_statement_correct = statement_correct or some_statement_correct

    assert(some_statement_correct)

def test_process_increase_amount():
    #Test extraction of IncreaseAmount statements from a text description
    #of a substance increasing the expression of some gene.

    s = 'BRAF increases the expression of p53.'

    tp = tees_process_text(s)
    statements = tp.statements

    # Exactly one statement should have been extracted from the provided text
    assert(len(statements) == 1)
    statement0 = statements[0]

    # The statement should be of type IncreaseAmount
    assert(isinstance(statement0, IncreaseAmount))

    # The subject and object in the statement should correspond with what is
    # in the text
    subj0 = statement0.subj.db_refs['TEXT']
    assert(subj0 == 'BRAF')
    #
    obj0 = statement0.obj.db_refs['TEXT']
    assert(obj0 == 'p53')

def test_process_decrease_amount():
    #Test extraction of DecreaseAmount statements from a text description
    #of a substance decreasing the expression of some gene.

    s = 'BRAF decreases the expression of p53.'
    tp = tees_process_text(s)
    statements = tp.statements

    # Exactly one statement should have been etracted from the provided text
    assert(len(statements) == 1)
    statement0 = statements[0]

    # The statement should be of type DecreaseAmount
    assert(isinstance(statement0, DecreaseAmount))

    # The subject and object in the statement should correspond with what is
    # in the text
    subj0 = statement0.subj.db_refs['TEXT']
    assert(subj0 == 'BRAF')
    #
    obj0 = statement0.obj.db_refs['TEXT']
    assert(obj0 == 'p53')

def test_process_bind():
    #Test extracting of Complex statement from a text description of
    #substances binding to each other.

    s = 'BRAF binds to p53.'
    tp = tees_process_text(s)
    statements = tp.statements

    # There should be exactly one statement extracted
    assert(len(statements) == 1)
    statement0 = statements[0]

    # The statement should be of type Complex
    assert(isinstance(statement0, Complex))

    # The members of the complex should correspond to the two substances
    # mentioned in the text
    members = statement0.members
    members = [m.db_refs['TEXT'] for m in members]
    assert(len(members) == 2)
    assert('BRAF' in members)
    assert('p53' in members)

