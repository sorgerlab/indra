Using natural language to build models
======================================

In this tutorial we build a simple model using natural language,
and export it into different formats.


Read INDRA Statements from a natural language string
----------------------------------------------------

First we import INDRA's API to the TRIPS reading system. We then define a block
of text which serves as the description of the mechanism to be modeled in the
`model_text` variable. Finally, `indra.sources.trips.process_text` is called
which sends a request to the TRIPS web service, gets a response and processes
the extraction knowledge base to obtain a list of INDRA Statements

.. Update code in tests/test_docs_code.py as well

.. ipython:: python

    from indra.sources import trips

    model_text = 'MAP2K1 phosphorylates MAPK1 and DUSP6 dephosphorylates MAPK1.'
    tp = trips.process_text(model_text)

At this point `tp.statements` should contain 2 INDRA Statements:
a Phosphorylation Statement and a Dephosphorylation Statement. Note that the
evidence sentence for each Statement is propagated:

.. Update code in tests/test_docs_code.py as well

.. ipython:: python

    for st in tp.statements:
        print('%s with evidence "%s"' % (st, st.evidence[0].text))

Assemble the INDRA Statements into a rule-based executable model
----------------------------------------------------------------

We next use INDRA's PySB Assembler to automatically assemble a rule-based model
representing the biochemical mechanisms described in `model_text`. First a
PysbAssembler object is instantiated, then the list of INDRA Statements is
added to the assembler. Finally, the assembler's `make_model` method is called
which assembles the model and returns it, while also storing it in `pa.model`.
Notice that we are using `policies='two_step'` as an argument of `make_model`.
This directs the assemble to use rules in which enzymatic catalysis is modeled
as a two-step process in which enzyme and substrate first reversibly bind and
the enzyme-substrate complex produces and releases a product irreversibly.

.. Update code in tests/test_docs_code.py as well

.. ipython:: python

    from indra.assemblers.pysb import PysbAssembler

    pa = PysbAssembler()
    pa.add_statements(tp.statements)
    pa.make_model(policies='two_step')

At this point `pa.model` contains a PySB model object with 3 monomers,

.. Update code in tests/test_docs_code.py as well

.. ipython:: python

    for monomer in pa.model.monomers:
        print(monomer)

6 rules,

.. Update code in tests/test_docs_code.py as well

.. ipython:: python

    for rule in pa.model.rules:
        print(rule)

and 9 parameters (6 kinetic rate constants and 3 total protein amounts) that
are set to nominal but plausible values,

.. Update code in tests/test_docs_code.py as well

.. ipython:: python

    for parameter in pa.model.parameters:
        print(parameter)

The model also contains extensive annotations that tie the monomers to database
identifiers and also annotate the semantics of each component of each rule.

.. Update code in tests/test_docs_code.py as well

.. ipython:: python

    for annotation in pa.model.annotations:
        print(annotation)

..  Set the model to a particular cell line context
    -----------------------------------------------
    We can use INDRA's contextualization module which is built into the
    PysbAssembler to set the amounts of proteins in the model to total amounts
    measured (or estimated) in a given cancer cell line. In this example,
    we will use the  `A375` melanoma cell line to set the total amounts of
    proteins in the model.
.. Update code in tests/test_docs_code.py as well
    .. ipython:: python
        pa.set_context('A375_SKIN')
    At this point the PySB model has total protein amounts set consistent with the
    A375 cell line:
.. Update code in tests/test_docs_code.py as well
    .. ipython:: python
        for monomer_pattern, parameter in pa.model.initial_conditions:
            print('%s = %d' % (monomer_pattern, parameter.value))

Exporting the model into other common formats
---------------------------------------------
From the assembled PySB format it is possible to export the model into other
common formats such as SBML, BNGL and Kappa. One can also generate a Matlab or
Mathematica script with ODEs corresponding to the model.

.. Update code in tests/test_docs_code.py as well

::

    pa.export_model('sbml')
    pa.export_model('bngl')

One can also pass a file name argument to the `export_model` function to save
the exported model directly into a file:

.. Update code in tests/test_docs_code.py as well

::

    pa.export_model('sbml', 'example_model.sbml')
