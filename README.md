[![Build Status](https://travis-ci.org/sorgerlab/indra.svg?branch=travis_ci)](https://travis-ci.org/sorgerlab/indra) [![Documentation Status](https://readthedocs.org/projects/indra/badge/?version=latest)](https://indra.readthedocs.io/en/latest/?badge=latest)

# INDRA

<img align="left" src="/doc/indra_logo.png?raw=True" width="300" height="224" />

INDRA (Integrated Network and Dynamical Reasoning Assembler) is an automated
model assembly system, originally developed for molecular systems biology and
currently being generalized to other domains. INDRA draws on natural language
processing systems and structured databases to collect mechanistic and causal
assertions, represents them in a standardized form (INDRA Statements), and
assembles them into various modeling formalisms including causal graphs and
dynamical models. INDRA also provides knowledge assembly procedures that
operate on INDRA Statements and correct certain errors, find and resolve
redundancies, infer missing information, filter to a scope of interest and
assess belief.

### Knowledge sources

INDRA is currently integrated with the following natural language processing
systems:
- [TRIPS/DRUM](http://trips.ihmc.us/parser/cgi/drum) - for biology
- [REACH](https://github.com/clulab/reach) - for biology
- [Sparser](https://github.com/ddmcdonald/sparser) - for biology
- [TEES](https://github.com/jbjorne/TEES) - for biology
- [Eidos](https://github.com/clulab/eidos) - general purpose

and can collect information from these databases:
- [Pathway Commons database](http://pathwaycommons.org/) or any source
    using the [BioPAX](http://www.biopax.org/) format
- [BEL Large Corpus](https://github.com/OpenBEL/) or any source using the
    [BEL](https://github.com/OpenBEL/) format
- [SIGNOR](https://signor.uniroma2.it/)

These input modules (available in `indra.sources`) all produce INDRA
Statements.

### Output model assemblers

INDRA also provides several model output assemblers that take INDRA Statements
as input. INDRA can assemble into the following modeling formalisms
- Detailed mechanistic, executable models in [PySB](http://pysb.org/)
    which can further be exported into SBML, SBGN, Kappa, and BNGL.
- Directed causal networks in
    [SIF](http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats), 
    [NDEx/CX](http://www.home.ndexbio.org/data-model/), 
    [Cytoscape.js](http://js.cytoscape.org/), and
    [Graphviz](https://www.graphviz.org/) formats.
- English language (a human-readable summary of the information
    collected and assembled by INDRA)

### Internal knowledge assembly

The internal assembly steps of INDRA are exposed in the
[indra.tools.assemble_corpus](http://indra.readthedocs.io/en/latest/modules/tools/index.html#module-indra.tools.assemble_corpus) 
submodule. This submodule contains functions that
take Statements as input and produce processed Statements as output. They can
be composed to form an assembly pipeline connecting knowledge collected from
sources with an output model.

INDRA also contains utility modules to access literature content (e.g. PubMed),
ontological information (e.g. UniProt, HGNC), and other resources.

## Citation

[From word models to executable models of signal transduction
using automated assembly](http://msb.embopress.org/content/13/11/954),
Molecular Systems Biology (2017)

## Documentation

Documentation is available at
[http://indra.readthedocs.io](http://indra.readthedocs.io).


## Installation

For detailed installation instructions,
[see the documentation](http://indra.readthedocs.io/en/latest/installation.html).

INDRA works with both Python 2 and 3 (tested with 2.7 and 3.5).

The preferred way to install INDRA is by pointing pip to the source repository
as

    $ pip install git+https://github.com/sorgerlab/indra.git

or by cloning the repository and then using pip to install the package as

    $ git clone https://github.com/sorgerlab/indra.git
    $ cd indra
    $ pip install .

You can also install INDRA by cloning this repository and running setup.py
as

    $ git clone https://github.com/sorgerlab/indra.git
    $ cd indra
    $ python setup.py install

Releases of INDRA are also available on
[PyPI](https://pip.pypa.io/en/latest/installing/), you can install the latest
release as

    $ pip install indra

However, releases will usually be behind the latest code available in this
repository.

INDRA depends on a few standard Python packages. These packages are installed by
either setup method (using pip or running setup.py install).
For certain modules and use cases, other dependencies may be needed,
which are described in detail in the
[documentation](http://indra.readthedocs.io/en/latest/installation.html).

## Using INDRA

In this example INDRA assembles a PySB model from the natural language
description of a mechanism via the [TRIPS reading web
service](http://trips.ihmc.us/parser/cgi/drum).

```python
from indra.sources import trips
from indra.assemblers import PysbAssembler
pa = PysbAssembler()
# Process a natural language description of a mechanism
trips_processor = trips.process_text('MEK2 phosphorylates ERK1 at Thr-202 and Tyr-204')
# Collect extracted mechanisms in PysbAssembler
pa.add_statements(trips_processor.statements)
# Assemble the model
model = pa.make_model(policies='two_step')
```

INDRA also provides an interface for the
[REACH](http://agathon.sista.arizona.edu:8080/odinweb/) natural language
processor. In this example, a full paper from [PubMed
Central](http://www.ncbi.nlm.nih.gov/pmc/) is processed. The paper's PMC ID is
[PMC3717945](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3717945/).

```python
from indra.sources import reach
# Process the neighborhood of BRAF and MAP2K1
reach_processor = reach.process_pmc('3717945')
# At this point, reach_processor.statements contains a list of INDRA statements
# extracted from the PMC paper.
```

Next we look at an example of reading the 10 most recent PubMed abstracts on
BRAF and collecting the results in INDRA statements.

```python
from indra.sources import reach
from indra.literature import pubmed_client
# Search for 10 most recent abstracts in PubMed on 'BRAF'
pmids = pubmed_client.get_ids('BRAF', retmax=10)
all_statements = []
for pmid in pmids:
    abs = pubmed_client.get_abstract(pmid)
    if abs is not None:
        reach_processor = reach.process_text(abs)
        if reach_processor is not None:
            all_statements += reach_processor.statements
# At this point, the all_statements list contains all the statements
# extracted from the 10 abstracts.
```

The next example shows querying the [BEL large
corpus](http://public.ndexbio.org/#/network/9ea3c170-01ad-11e5-ac0f-000c29cb28fb)
network through [NDEx](http://ndexbio.org) for a neighborhood of a given list
of proteins using their HGNC gene names.

```python
from indra.sources import bel
# Process the neighborhood of BRAF and MAP2K1
bel_processor = bel.process_ndex_neighborhood(['BRAF', 'MAP2K1'])
# At this point, bel_processor.statements contains a list of INDRA statements
# extracted from the neihborhood query.
```

Next, we look at an example of querying the [Pathway Commons
database](http://pathwaycommons.org) for paths between two lists of proteins.
Note: see installation notes above for installing jnius, which is required for
using the BioPAX API of INDRA.

```python
from indra.sources import biopax
# Process the neighborhood of BRAF and MAP2K1
biopax_processor = biopax.process_pc_pathsfromto(['BRAF', 'RAF1'], ['MAP2K1', 'MAP2K2'])
# At this point, biopax_processor.statements contains a list of INDRA 
# Statements extracted from the paths-from-to query.
```

