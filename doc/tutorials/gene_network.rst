Assembling everything known about a particular gene
===================================================

Assume you are interested in collecting all mechanisms that a particular gene
is involved in. Using INDRA, it is possible to collect everything curated
about the gene in pathway databases and then read all the accessible literature
discussing the gene of interest. This knowledge is aggregated as a set of
INDRA Statements which can then be assembled into several different model
and network formats and possibly shared online.

For the sake of this example, assume that the gene of interest is H2AX.

It is important to use the standard HGNC gene symbol of the gene throughout the
example (this information is available on http://www.genenames.org/ or
http://www.uniprot.org/) - abritrary synonyms will not work!

Collect mechanisms from PathwayCommons and the BEL Large Corpus
---------------------------------------------------------------

We first collect Statements from the PathwayCommons database via INDRA's
BioPAX API and then collect Statements from the BEL Large Corpus via INDRA's
BEL API.

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code-block:: python

    from indra.tools.gene_network import GeneNetwork

    gn = GeneNetwork(['H2AX'])
    biopax_stmts = gn.get_biopax_stmts()
    bel_stmts = gn.get_bel_stmts()

at this point `biopax_stmts` and `bel_stmts` are two lists of INDRA Statements.

Collect a list of publications that discuss the gene of interest
----------------------------------------------------------------

We next use INDRA's literature client to find PubMed IDs (PMIDs) that discuss
the gene of interest. To find articles that are annotated with the given gene,
INDRA first looks up the Entrez ID corresponding to the gene name and then
finds associated publications.

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code-block:: python

    from indra import literature

    pmids = literature.pubmed_client.get_ids_for_gene('H2AX')

The variable `pmids` now contains a list of PMIDs associated with the gene.

Get the abstracts corresponding to the publications
---------------------------------------------------

Next we use INDRA's literature client to fetch the abstracts corresponding to
the PMIDs we have just collected. The client also returns other content
types, like xml, for full text (if available). Here we cut the list of PMIDs
short to just the first 10 IDs that contain abstracts to make the processing
faster.

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code-block:: python

    from indra import literature

    paper_contents = {}
    for pmid in pmids:
        content, content_type = literature.get_full_text(pmid, 'pmid')
        if content_type == 'abstract':
            paper_contents[pmid] = content
        if len(paper_contents) == 10:
            break

We now have a dictionary called `paper_contents` which stores the content for
each PMID we looked up. While the abstracts are in plain text format,
some content is sometimes returned in different either PMC NXML or Elsevier
XML format. To process
XML from different sources, some example are:
`INDRA Reach API <../modules/sources/reach/index.html#indra.sources.reach
.api.process_nxml_str>`_ or the
`INDRA Elsevier client <../modules/literature/index.html#module-indra
.literature.elsevier_client>`_.

Read the content of the publications
------------------------------------

We next run the REACH reading system on the publications. Here we assume
 that the REACH web service is running locally and is available at
 `http://localhost:8080` (the default web service endpoints for processing
 text and nxml are available as importable variables e.g., `local_text_url`.
 To get started wtih this, see method 1 listed in <`INDRA Reach API
 <../modules/sources/reach/index.html#indra.sources.reach>`_ documentation.

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code-block:: python

    from indra.sources import reach

    literature_stmts = []
    for pmid, content in paper_contents.items():
        rp = reach.process_text(content, url=reach.local_text_url)
        literature_stmts += rp.statements
    print('Got %d statements' % len(literature_stmts))

The list `literature_stmts` now contains the results of all the statements
that were read.

Combine all statements and run pre-assembly
-------------------------------------------

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code-block:: python

    from indra.tools import assemble_corpus as ac

    stmts = biopax_stmts + bel_stmts + literature_stmts

    stmts = ac.map_grounding(stmts)
    stmts = ac.map_sequence(stmts)
    stmts = ac.run_preassembly(stmts)

At this point `stmts` contains a list of Statements with
`grounding <../modules/preassembler/grounding_mapper.html>`_, having been
mapped according to INDRA's built in grounding map and disambiguation
features, amino acid sites having been
`mapped <../modules/preassembler/site_mapper.html>`_,
duplicates combined, and hierarchically subsumed variants of statements hidden.
It is possible to run other assembly steps and filters on the results such as
to keep only human genes, remove Statements with ungrounded genes, or to
keep only certain types of interactions. You can find more assembly steps that
can be included in your pipeline in the `Assemble Corpus documentation
<../modules/tools/index.html#module-indra.tools.assemble_corpus>`_.
You can also read more about the pre-assembly
process in the
`preassembly module documentation <../modules/preassembler/index.html>`_ and
in the `GitHub documentation
<https://github.com/sorgerlab/indra#internal-knowledge-assembly>`_

Assemble the statements into a network model
--------------------------------------------

CX Network Model
~~~~~~~~~~~~~~~~

We can assemble the statements into e.g., a CX network model:

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code-block:: python

    from indra.assemblers.cx import CxAssembler
    from indra.databases import ndex_client

    cxa = CxAssembler(stmts)
    cx_str = cxa.make_model()

We can now upload this network to the Network Data Exchange (NDEx).

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code-block:: python

    ndex_cred = {'user': 'myusername', 'password': 'xxx'}
    network_id = ndex_client.create_network(cx_str, ndex_cred)
    print(network_id)

IndraNet Model
~~~~~~~~~~~~~~

Another network model that can assembled is the IndraNet graph which is a
light-weight networkx derived object.

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code:: python

    from indra.assemblers.indranet import IndraNetAssembler
    indranet_assembler = IndraNetAssembler(statements=stmts)
    indranet = indranet_assembler.make_model()

Since the IndraNet class is a child class of a networkx Graph, one can use
networkx's algorithms:

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code:: python

    import networkx as nx
    paths = nx.single_source_shortest_path(G=indranet, source='H2AX',
                                           cutoff=1)

Executable PySB Model
~~~~~~~~~~~~~~~~~~~~~

An executable PySB model can be assembled with the PySB assembler:

.. Update code in tests/test_docs_code.py:test_gene_network as well

.. code:: python

    from indra.assemblers.pysb import PysbAssembler
    pysb = PysbAssembler(statements=stmts)
    pysb_model = pysb.make_model()

Read more about PySB models in the `PySB documentation <http://pysb.org/>`_
and look into the `natural language modeling tutorial <nl_modeling.html>`_
which uses PySB models.

Read more about all assembly output formats in the `README <https://github
.com/sorgerlab/indra#output-model-assemblers>`_ and in the `module
references <https://indra.readthedocs.io/en/latest/modules/assemblers/index
.html>`_.
