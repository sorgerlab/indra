from nose.plugins.attrib import attr
from indra.statements import *
from indra.sources.tees import api
from indra.sources.tees.processor import TEESProcessor

_multiprocess_can_split_ = False
_multiprocess_shared_ = False


@attr('slow')
def test_process_phosphorylation():
    # Test the extraction of phosphorylation with a simple example.
    s = 'Ras leads to the phosphorylation of Braf.'
    tp = api.process_text(s)
    statements = tp.statements

    # Should only extract one statement
    assert len(statements) == 1

    # The statement should be a Phosphorylation statement
    statement = statements[0]
    assert isinstance(statement, Phosphorylation)

    # Check that the enzyme and substrate of the statement match the text
    enz = statement.enz
    assert enz.db_refs['TEXT'] == 'Ras'
    #
    sub = statement.sub
    assert sub.db_refs['TEXT'] == 'Braf'

    # There should be an evidence object (properties of evidence tested in
    # other tests)
    assert len(statements[0].evidence) == 1


@attr('slow')
def test_process_dephosphorylation():
    # Test the extraction of a dephosphorylation sentence. This sentence is
    # processed into two INDRA statements, at least one of which is correct.

    # Text statement describing phosphorylation
    s = 'Here we show that febuxostat suppresses LPS-induced MCP-1 ' + \
            'production and mRNA expression via activating ' + \
            'MAPK phosphatase-1 (MKP-1) which, in turn, leads to ' + \
            'dephosphorylation and inactivation of JNK in macrophages.'
    tp = api.process_text(s)
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

        # There should be an evidence object (properties of evidence tested in
        # other tests)
        assert len(statements[0].evidence) == 1

        # If yes to both, the correct statement was among those extracted
        statement_correct = type_correct and entities_correct
        some_statement_correct = statement_correct or some_statement_correct

    assert some_statement_correct


@attr('slow')
def test_process_increase_amount():
    # Test extraction of IncreaseAmount statements from a text description
    # of a substance increasing the expression of some gene.

    s = 'BRAF increases the expression of p53.'

    tp = api.process_text(s)
    statements = tp.statements

    # Exactly one statement should have been extracted from the provided text
    assert len(statements) == 1
    statement0 = statements[0]

    # The statement should be of type IncreaseAmount
    assert isinstance(statement0, IncreaseAmount)

    # The subject and object in the statement should correspond with what is
    # in the text
    subj0 = statement0.subj.db_refs['TEXT']
    assert subj0 == 'BRAF'
    #
    obj0 = statement0.obj.db_refs['TEXT']
    assert obj0 == 'p53'

    # There should be an evidence object (properties of evidence tested in
    # other tests)
    assert len(statements[0].evidence) == 1


@attr('slow')
def test_process_decrease_amount():
    # Test extraction of DecreaseAmount statements from a text description
    # of a substance decreasing the expression of some gene.

    s = 'BRAF decreases the expression of p53.'
    tp = api.process_text(s)
    statements = tp.statements

    # Exactly one statement should have been extracted from the provided text
    assert len(statements) == 1
    statement0 = statements[0]

    # The statement should be of type DecreaseAmount
    assert isinstance(statement0, DecreaseAmount)

    # The subject and object in the statement should correspond with what is
    # in the text
    subj0 = statement0.subj.db_refs['TEXT']
    assert subj0 == 'BRAF'
    #
    obj0 = statement0.obj.db_refs['TEXT']
    assert obj0 == 'p53'

    # There should be an evidence object (properties of evidence tested in
    # other tests)
    assert len(statements[0].evidence) == 1


@attr('slow')
def test_process_bind():
    # Test extracting of Complex statement from a text description of
    # substances binding to each other.

    s = 'BRAF binds to p53.'
    tp = api.process_text(s)
    statements = tp.statements

    # There should be exactly one statement extracted
    assert len(statements) == 1
    statement0 = statements[0]

    # The statement should be of type Complex
    assert isinstance(statement0, Complex)

    # The members of the complex should correspond to the two substances
    # mentioned in the text
    members = statement0.members
    members = [m.db_refs['TEXT'] for m in members]
    assert len(members) == 2
    assert 'BRAF' in members
    assert 'p53' in members

    # is_direct should be true in evidence for bind statement
    assert len(statement0.evidence) == 1
    assert statement0.evidence[0].epistemics['direct']


@attr('slow')
def test_evidence_text():
    # Test the ability of the processor to extract which sentence in particular
    # lead to the creation of the INDRA statement, amongst a corpus of text
    # not related to the statement's mechanism.

    # Corpus containing exactly one biological mechanism
    corpus = """Why did the cows return to the marijuana field? It was the pot
    calling the cattle back. Why do cows have hooves instead of feet? Because
    they lactose. When making non-dairy butter, there is little margarine for
    error. Ras leads to the phosphorylation of Raf. Do ghost cows say "moo" or
    "boo"? The surprising fact is they say "moo", but with a rising vibrato
    tone instead of an elongated one."""

    # Process the corpus
    tp = api.process_text(corpus)
    statements = tp.statements

    # Only one of the sentences was related to a biological mechanism
    assert len(statements) == 1
    statement0 = statements[0]

    # Make sure it got the right sentence
    text = statement0.evidence[0].text
    print('Statement text:"', text + '"')
    assert text == 'Ras leads to the phosphorylation of Raf.'


#@attr('slow')
def test_evidence_pmid():
    # Test whether the pmid provided to the TEES processor is put into the
    # statement's evidence

    pmid = '42'

    # Process some text
    s = 'BRAF binds to p53.'
    tp = api.process_text(s, pmid)
    statements = tp.statements

    # Verify that the pmid was put in the right place
    assert len(statements[0].evidence) == 1
    assert statements[0].evidence[0].pmid == pmid


def test_ground():
    mek = Agent('Mek', db_refs={'TEXT': 'MEK'})
    erk = Agent('Erk1', db_refs={'TEXT': 'Erk1'})
    stmt = Phosphorylation(mek, erk)
    TEESProcessor.ground_statements([stmt])
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.enz.db_refs['FPLX'] == 'MEK'
    assert stmt.sub.name == 'MAPK3'
    assert stmt.sub.db_refs['HGNC'] == '6877'
    assert stmt.sub.db_refs['UP'] == 'P27361'
