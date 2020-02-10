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

    from indra.sources import bel

Similarly, each model output assembler has its own submodule under
`indra.assemblers` with the assembler class accessible at the submodule
level, so they can be imported as, for instance,

.. code:: python

    from indra.assemblers.pysb import PysbAssembler

To get a detailed overview of INDRA's submodule structure, take a look at
the :ref:`indra_modules_ref`.

Basic usage examples
--------------------

Here we show some basic usage examples of the submodules of INDRA. More complex
usage examples are shown in the Tutorials section.

Reading a sentence with TRIPS
`````````````````````````````
In this example, we read a sentence via INDRA's TRIPS submodule to produce
an INDRA Statement.

.. code:: python

    from indra.sources import trips
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

    from indra.sources import reach
    reach_processor = reach.process_pmc('3717945')

The `reach_processor` object has a `statements` attribute which contains a list
of INDRA Statements extracted from the paper.

Getting the neighborhood of proteins from the BEL Large Corpus
``````````````````````````````````````````````````````````````
In this example, we search the neighborhood of the KRAS and BRAF proteins in
the BEL Large Corpus.

.. code:: python

    from indra.sources import bel
    bel_processor = bel.process_pybel_neighborhood(['KRAS', 'BRAF'])

The `bel_processor` object has a `statements` attribute which contains a list
of INDRA Statements extracted from the queried neighborhood.

Getting paths between two proteins from PathwayCommons (BioPAX)
```````````````````````````````````````````````````````````````
In this example, we search for paths between the BRAF and MAPK3 proteins in the
PathwayCommons databases using INDRA's BioPAX API. Note that this example will
only work if all dependencies of the indra.sources.biopax module are
installed.

See the `Installation instructions <installation.html>`_ for more details.

.. code:: python

    from indra.sources import biopax
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

Assembling a PySB model and exporting to SBML
`````````````````````````````````````````````
In this example, assume that we have already collected a list of INDRA Statements
from any of the input sources and that this list is called `stmts`. We will
instantiate a PysbAssembler, which produces a PySB model from INDRA Statements.

.. code:: python

    from indra.assemblers.pysb import PysbAssembler
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()

Here the `model` variable is a PySB Model object representing a rule-based
executable model, which can be further manipulated, simulated, saved and exported
to other formats.

For instance, exporting the model to `SBML <http://sbml.org>`_ format can
be done as

.. code:: python

    sbml_model = pa.export_model('sbml')

which gives an SBML model string in the `sbml_model` variable, or as

.. code:: python

    pa.export_model('sbml', file_name='model.sbml')

which writes the SBML model into the `model.sbml` file. Other formats for export
that are supported include BNGL, Kappa and Matlab. For a full list, see the
`PySB export module
<http://docs.pysb.org/en/latest/modules/export/index.html>`_.

Exporting Statements As Graph
`````````````````````````````
In this example we again assume that there already exists a variable called
`stmts`, containing a list of statements. We will import the
`IndraNetAssembler` that produces a multidigraph representations of a list
of statements.

.. code:: python

    from indra.assemblers.indranet import IndraNetAssembler
    indranet_assembler = IndraNetAssembler(statements=stmts)
    indranet = indranet_assembler.make_model()

The indranet object is an instance of a childclass of a Networkx graph,
making all networkx graph methods available for the indranet object. Each
edge in the graph represents one statement and has an edge dictionary with
meta data about the statement.

The indranet graph has methods to map it to other graph types. Here we
export it to a signed graph:

.. code:: python

    signed_graph = indranet.to_signed_graph()

Read more about the `IndraNetAssembler` in the `documentation
<modules/assemblers/indranet_assembler.html>`_.
