Getting started with INDRA
==========================

Importing INDRA and its modules
-------------------------------
INDRA can be imported and used in a Python script or interactively
in a Python shell. Note that similar to some other packages (e.g scipy), INDRA
doesn't automatically import all its submodules, so
`import indra` is not enough to access its submodules.
Rather, one has to explicitly import each submodule that is needed.
For example to access the BEL API, one has to

.. code:: python

    from indra import bel

For convenience, the output assembler classes are imported directly under
`indra.assemblers` so they can be imported as, for instance,

.. code:: python

    from indra.assemblers import PysbAssembler

To get a detailed overview of INDRA's submodule structure, take a look at the :ref:`indra_modules_ref`.

Basic usage examples
--------------------

Here we show some basic usage examples of the submodules of INDRA. More complex
usage examples are shown in the Tutorials section.

Reading a sentence with TRIPS
`````````````````````````````
In this example, we read a sentence via INDRA's TRIPS submodule to produce
an INDRA Statement.

.. code:: python

    from indra import trips
    sentence = 'MAP2K1 phosphorylates MAPK3 at Thr-202 and Tyr-204'
    trips_processor = trips.process_text(sentence)

The `trips_processor` object has a `statements` attribute which contains a list
of INDRA Statements extracted from the sentence.

Reading a PubMed Central article with REACH
```````````````````````````````````````````
In this example, a full paper from `PubMed
Central <http://www.ncbi.nlm.nih.gov/pmc/>`_ is processed. The paper's PMC ID is
`PMC3717945 <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3717945/>`_.

.. code:: python

    from indra import reach
    reach_processor = reach.process_pmc('3717945')

The `reach_processor` object has a `statements` attribute which contains a list
of INDRA Statements extracted from the paper.

Getting the neighborhood of proteins from the BEL Large Corpus
``````````````````````````````````````````````````````````````
In this example, we search the neighborhood of the KRAS and BRAF proteins in
the BEL Large Corpus.

.. code:: python

    from indra import bel
    bel_processor = bel.process_ndex_neighborhood(['KRAS', 'BRAF'])

The `bel_processor` object has a `statements` attribute which contains a list
of INDRA Statements extracted from the queried neighborhood.

Getting paths between two proteins from PathwayCommons (BioPAX)
```````````````````````````````````````````````````````````````
In this example, we search for paths between the BRAF and MAPK3 proteins in
the PathwayCommons databases using INDRA's BioPAX API. Note that this example
will only work if all dependencies of the indra.biopax module are installed.
See the `Installation instructions <installation.html>`_ for more details.

.. code:: python

    from indra import biopax
    proteins = ['BRAF', 'MAPK3']
    limit = 2
    biopax_processor = biopax.process_pc_pathsbetween(proteins, limit)

We passed the second argument `limit = 2`, which defines the upper limit on
the length of the paths that are searched. By default the limit is 1.
The `biopax_processor` object has a `statements` attribute which contains a list
of INDRA Statements extracted from the queried paths.

Constructing INDRA Statements manually
``````````````````````````````````````
It is possible to construct INDRA Statements manually or in scripts. The following
is a basic example in which we instantiate a Phosphorylation Statement between
BRAF and MAP2K1.

.. code:: python

    from indra.statements import Phosphorylation, Agent
    braf = Agent('BRAF')
    map2k1 = Agent('MAP2K1')
    stmt = Phosphorylation(braf, map2k1)

Assembling a PySB model
```````````````````````
In this example, assume that we have already collected a list of INDRA Statements
from any of the input sources and that this list is called `stmts`. We will
instantiate a PysbAssembler, which produces a PySB model from INDRA Statements.

.. code:: python

    from indra.assemblers import PysbAssembler
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()

Here the `model` variable is a PySB Model object representing a rule-based
executable model, which can be further manipulated, simulated, saved and exported
to other formats.
