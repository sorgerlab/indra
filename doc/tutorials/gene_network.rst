Assembling everything known about a particular gene
===================================================

Assume you are interested in collecting all mechanisms that a particular gene
is involved in. Using INDRA, it is possible to collect everything curated
about the gene in pathway databases and then read all the accessible literature
discussing the gene of interest. This knowledge is aggregated as a set of
INDRA Statements which can then be assembled into several different model
and network formats and possibly shared online.

For the sake of example, assume that the gene of interest is STING1.

It is important to use the standard HGNC gene symbol of the gene throughout the
example (this information is available on http://www.genenames.org/ or
http://www.uniprot.org/) - abritrary synonyms will not work!

Collect mechanisms from PathwayCommons and the BEL Large Corpus
---------------------------------------------------------------

We first collect Statements from the PathwayCommons database via INDRA's
BioPAX API and then collect Statements from the BEL Large Corpus via INDRA's
BEL API.

.. code-block:: python

    from indra.tools.gene_network import GeneNetwork

    gn = GeneNetwork(['STING1'])
    biopax_stmts = gn.get_biopax_stmts()
    bel_stmts = gn.get_bel_stmts()

at this point `biopax_stmts` and `bel_stmts` are two lists of INDRA Statements.

Collect a list of publications that discuss the gene of interest
----------------------------------------------------------------

We next use INDRA's literature client to find PubMed IDs (PMIDs) that discuss
the gene of interest. To find articles that are annotated with the given gene,
INDRA first looks up the Entrez ID corresponding to the gene name and then
finds associated publications.

.. code-block:: python

    from indra import literature

    pmids = literature.pubmed_client.get_ids_for_gene('STING1')

The variable `pmid` now contains a list of PMIDs associated with the gene.

Get the full text or abstract corresponding to the publications
---------------------------------------------------------------

Next we use INDRA's literature client to fetch the full text (if available) or
the abstract corresponding to the PMIDs we have just collected.

.. code-block:: python

    from indra import literature

    paper_contents = {}
    for pmid in pmids:
        content, content_type = literature.get_full_text(pmid, 'pmid')
        paper_contents[pmid] = (content, content_type)

We now have a dictionary called `paper_contents` which stores the content and
the content type of each PMID we looked up.

Read the content of the publications
------------------------------------

We next run the REACH reading system on the publications. Depending on the 
content type, different calls need to be made via INDRA's REACH API.

.. code-block:: python

    from indra import literature
    from indra.sources import reach

    read_offline = True

    literature_stmts = []
    for pmid, (content, content_type) in paper_contents.items():
        rp = None
        print('Reading %s' % pmid)
        if content_type == 'abstract':
            rp = reach.process_text(content, citation=pmid, offline=read_offline)
        elif content_type == 'pmc_oa_xml':
            rp = reach.process_nxml_str(content, offline=read_offline)
        elif content_type == 'elsevier_xml':
            txt = literature.elsevier_client.extract_text(content)
            if txt:
                rp = reach.process_text(txt, citation=pmid, offline=read_offline)
        if rp is not None:
            literature_stmts += rp.statements

The list `literature_stmts` now contains the results of all the statements
that were read.

Combine all statements and run pre-assembly
-------------------------------------------

.. code-block:: python

    from indra.tools import assemble_corpus

    stmts = biopax_stmts + bel_stmts + literature_stmts

    stmts = assemble_corpus.map_grounding(stmts)
    stmts = assemble_corpus.map_sequence(stmts)
    stmts = assemble_corpus.run_preassembly(stmts)

At this point `stmts` contains a list of Statements collected with grounding,
sequences having been mapped, duplicates combined and less specific variants
of statements hidden. It is possible to run other filters on the results such
as to keep only human genes, remove Statements with ungrounded genes, or
to keep only certain types of interactions.

Assemble the statements into a network model
--------------------------------------------

.. code-block:: python

    from indra.assemblers.cx import CxAssembler
    from indra.databases import ndex_client

    cxa = CxAssembler(stmts)
    cx_str = cxa.make_model()

We can now upload this network to the Network Data Exchange (NDEx).

.. code-block:: python

    ndex_cred = {'user': 'myusername', 'password': 'xxx'}
    network_id = ndex_client.create_network(cx_str, ndex_cred)
    print(network_id)
