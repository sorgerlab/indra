[![Build Status](https://travis-ci.org/sorgerlab/indra.svg?branch=travis_ci)](https://travis-ci.org/sorgerlab/indra) [![Documentation Status](https://readthedocs.org/projects/indra/badge/?version=latest)](https://indra.readthedocs.io/en/latest/?badge=latest)

# INDRA

<img align="left" src="/doc/indra_logo.png?raw=True" width="300" height="224" />

INDRA (Integrated Network and Dynamical Reasoning Assembler) is an automated
model assembly system, originally developed for molecular systems biology and
currently being generalized to other domains. INDRA draws on natural language
processing systems and structured databases to collect mechanistic and causal
assertions, represents them in a standardized form (INDRA Statements), and
assembles them into various modeling formalisms including causal graphs and
dynamical models.

Importantly, INDRA isn't merely an import/export hub or a multi-reader
integrator. Its core added value in
assembling knowledge sources into coherent models comes from
providing knowledge-level assembly procedures that
operate on INDRA Statements and correct certain errors, find and resolve
redundancies, infer missing information, filter to a scope of interest and
assess belief.

The detailed INDRA documentation is available at
[http://indra.readthedocs.io](http://indra.readthedocs.io).

### Knowledge sources

INDRA is currently integrated with the following natural language processing
systems and structured databases. These input modules (available in
`indra.sources`) all produce INDRA Statements.

General purpose causal relation reading systems:

| Reader     | Module                | Reference                                 |
|------------|-----------------------|-------------------------------------------|
| Eidos      | `indra.sources.eidos` | https://github.com/clulab/eidos           |
| TRIPS/CWMS | `indra.sources.cwms`  | http://trips.ihmc.us/parser/cgi/cwmsreader|
| Hume       | `indra.sources.hume`  | https://github.com/BBN-E/Hume             |
| Sofia      | `indra.sources.sofia` | https://sofia.worldmodelers.com/ui/       |

Biology-oriented reading systems:

| Reader     | Module                  | Reference                                       |
|------------|-------------------------|-------------------------------------------------|
| TRIPS/DRUM | `indra.sources.trips`   | http://trips.ihmc.us/parser/cgi/drum            |
| REACH      | `indra.sources.reach`   | https://github.com/clulab/reach                 |
| Sparser    | `indra.sources.sparser` | https://github.com/ddmcdonald/sparser           |
| TEES       | `indra.sources.tees`    | https://github.com/jbjorne/TEES                 |
| MedScan    | `indra.sources.medscan` | https://doi.org/10.1093/bioinformatics/btg207   |
| RLIMS-P    | `indra.sources.rlimsp`  | https://research.bioinformatics.udel.edu/rlimsp |
| ISI/AMR    | `indra.sources.isi`     | https://github.com/sgarg87/big_mech_isi_gg      |
| Geneways   | `indra.sources.geneways`| https://www.ncbi.nlm.nih.gov/pubmed/15016385    |

Biological pathway databases:

| Database / Exchange format | Module                     | Reference                                                       |
|----------------------------|----------------------------|-----------------------------------------------------------------|
| PathwayCommons / BioPax    | `indra.sources.biopax`     | http://pathwaycommons.org/ <br/> http://www.biopax.org/         |
| Large Corpus / BEL         | `indra.sources.bel`        | https://github.com/pybel/pybel <br/> https://github.com/OpenBEL |
| Signor                     | `indra.sources.signor`     | https://signor.uniroma2.it/                                     |
| BioGRID                    | `indra.sources.biogrid`    | https://thebiogrid.org/                                         |
| Target Affinity Spectrum   | `indra.sources.tas`        | https://doi.org/10.1101/358978                                  |
| LINCS small molecules      | `indra.sources.lincs_drug` | http://lincs.hms.harvard.edu/db/sm/                             |

Custom knowledge bases:

| Database / Exchange format | Module                        | Reference                            |
|----------------------------|-------------------------------|--------------------------------------|
| NDEx / CX                  | `indra.sources.ndex_cx`       | http://ndexbio.org                   |
| INDRA DB / INDRA Statements| `indra.sources.indra_db_rest` | https://github.com/indralab/indra_db |


### Output model assemblers

INDRA also provides several model output assemblers that take INDRA Statements
as input. The most sophisticated model assembler is the PySB assembler, which
implements a policy-guided automated assembly procedure of a rule-based
executable model (that can then be further compiled into other formats such as
SBML, Kappa, BNGL and SBGN to connect to a vast ecosystem of downstream tools)
from INDRA Statements. Several other model assembly
modules target various network formats for browsing, display, and
graph/structural analysis (CyJS, Graphviz, SBGN, CX, SIF). Finally, the English
assembler produces.

INDRA also supports extension by outside model assembly tools which take
INDRA Statements as input and produce models. One such example is Delphi
(https://github.com/ml4ai/delphi), which is a Dynamic Bayesian Network
model assembler. Similarly, outside tools that support INDRA Statements
can implement custom visualization methods, such as CauseMos, developed
by Uncharted Software (https://uncharted.software/).

Assemblers aimed at model-driven discovery and analysis:

| Modeling formalism / Exchange format           | Purpose                                              | Module                  | Reference           |
|------------------------------------------------|------------------------------------------------------|-------------------------|---------------------|
| PySB (-> SBML, SBGN, BNGL, Kappa, etc.)        | Detailed, mechanistic modeling, simulation, analysis | `indra.assemblers.pysb` | http://pysb.org     |
| PyBEL                                          | Causal analysis, visualization                       | `indra.assemblers.pybel`| https://github.com/pybel/pybel <br/> https://bel-commons.scai.fraunhofer.de/ |
| SIF                                            | Network analysis, logic modeling, visualization      | `indra.assemblers.sif`  | [SIF format](http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html#sif-format) |
| Figaro                                         | Bayesian network inference                           | `indra.assemblers.figaro` | https://github.com/p2t2/figaro/ |
| KAMI                                           | Knowledge aggregation of protein sites/states and Kappa modeling | `indra.assemblers.kami` | https://github.com/Kappa-Dev/KAMI |


Assemblers primarily aimed at visualization:

| Network / Exchange format                      | Purpose                                              | Module                  | Reference           |
|------------------------------------------------|------------------------------------------------------|-------------------------|---------------------|
| Causal Analysis Graph                          | General causal graph visualization                   | `indra.assemblers.cag`  |                     |
| CX                                             | Network browsing, versioning on NDEx                 | `indra.assemblers.cx`   | http://ndexbio.org  |
| Cytoscape JS                                   | Interactive Cytoscape JS network to embed in websites| `indra.assemblers.cyjs` | http://js.cytoscape.org/ |
| Graphviz                                       | Static PDF/PNG visualization with powerful automated layout using Graphviz | `indra.assemblers.graph` | https://www.graphviz.org/ |
| SBGN                                           | Visualization with Systems Biology Graphical Notation| `indra.assemblers.sbgn` | http://sbgn.org     |

Assemblers primarily aimed at expert curation and browsing:

| Output format                                  | Purpose                                                | Module                  | Reference           |
|------------------------------------------------|------------------------------------------------------  |-------------------------|---------------------|
| English language                               | Human-readable descriptions, reports, dialogue         | `indra.assemblers.english` |                  |
| HTML                                           | Web-based browsing, linking out to provenance, curation| `indra.assemblers.html` | [Curation tutorial](https://indra.readthedocs.io/<br/>en/latest/tutorials/html_curation.html) |
| TSV (Tab/Comma Separated Values)               | Spreadsheet-based browsing and curation                | `indra.assemblers.tsv`  |                     |
| Index Cards                                    | Custom JSON format for curating biological mechanisms  | `indra.assemblers.index_card` |               |

### Internal knowledge assembly

A key feature of INDRA is providing internal knowledge-assembly modules
that operate on INDRA Statements and perform the following tasks:
- Redundancy/subsumption/generalization/contradiction finding and resolution
with respect to an ontology with the Preassembler
(`indra.preassembler.Preassembler`)
- Belief calculation based on evidence using the BeliefEngine
(`indra.belief`)
- Mapping grounding between multiple ontologies
(`indra.preassembler.ont_mapper.OntMapper`)
- Grounding override and disambiguation
(`indra.preassembler.grounding_mapper.GroundingMapper`)
- Protein sequence mapping (`indra.preassembler.site_mapper.SiteMapper`)

The internal assembly steps of INDRA including the ones listed above, and also
a large collection of filters (filter by source, belief, the presence of
grounding information, semantic filters by entity role, etc.) are exposed
in the
[indra.tools.assemble_corpus](http://indra.readthedocs.io/en/latest/modules/tools/index.html#module-indra.tools.assemble_corpus) 
submodule. This submodule contains functions that
take Statements as input and produce processed Statements as output. They can
be composed to form an assembly pipeline connecting knowledge collected from
sources with an output model.

INDRA also contains utility modules to access literature content (e.g. PubMed),
ontological information (e.g. UniProt, HGNC), and other resources.

## Citation

Gyori B.M., Bachman J.A., Subramanian K., Muhlich J.L., Galescu L., Sorger P.K.
[From word models to executable models of signaling networks using automated
assembly](http://msb.embopress.org/content/13/11/954) (2017),
Molecular Systems Biology, 13, 954.

## Installation

For detailed installation instructions,
[see the documentation](http://indra.readthedocs.io/en/latest/installation.html).

INDRA currently supports Python 3.5+. The last release of INDRA compatible
with Python 2.7 was 1.10.

The preferred way to install INDRA is by pointing pip to the source repository
as

    $ pip install git+https://github.com/sorgerlab/indra.git

Releases of INDRA are also available on
[PyPI](https://pip.pypa.io/en/latest/installing/), you can install the latest
release as

    $ pip install indra

However, releases will usually be behind the latest code available in this
repository.

INDRA depends on a few standard Python packages. These packages are installed
by pip during setup.
For certain modules and use cases, other "extra" dependencies may be needed,
which are described in detail in the
[documentation](http://indra.readthedocs.io/en/latest/installation.html).


## INDRA REST API
A REST API for INDRA is available at http://api.indra.bio:8000 with
documentation at http://www.indra.bio/rest_api/docs. Note that the REST API
is ideal for prototyping and for building light-weight web apps, but should
not be used for large reading and assembly workflows.


## INDRA Docker
INDRA is available as a Docker image on Dockerhub and can be pulled as

```
docker pull labsyspharm/indra
```

You can run the INDRA REST API using the container as
```
docker run -id -p 8080:8080 --entrypoint python labsyspharm/indra /sw/indra/rest_api/api.py
```

To build the image locally, there are currently two Dockerfiles for
INDRA and its dependencies. They are available in the following repositories:
- https://github.com/indralab/indra_docker
- https://github.com/indralab/indra_deps_docker

## Using INDRA

In this example INDRA assembles a PySB model from the natural language
description of a mechanism via the [TRIPS reading web
service](http://trips.ihmc.us/parser/cgi/drum).

```python
from indra.sources import trips
from indra.assemblers.pysb import PysbAssembler
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
```
At this point, `reach_processor.statements` contains a list of INDRA statements
extracted from the PMC paper.

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
```
At this point, the `all_statements` list contains all the statements
extracted from the 10 abstracts.

The next example shows querying the [BEL large
corpus](http://public.ndexbio.org/#/network/9ea3c170-01ad-11e5-ac0f-000c29cb28fb)
network for a neighborhood of a given list of proteins using their
HGNC gene names. This example performs the query via PyBEL.

```python
from indra.sources import bel
# Process the neighborhood of BRAF and MAP2K1
bel_processor = bel.process_pybel_neighborhood(['BRAF', 'MAP2K1'])
```
At this point, `bel_processor.statements` contains a list of INDRA statements
extracted from the neighborhood query.

Next, we look at an example of querying the [Pathway Commons
database](http://pathwaycommons.org) for paths between two lists of proteins.
Note: see installation notes above for installing pyjnius, which is required
for using the BioPAX API of INDRA.

```python
from indra.sources import biopax
# Process the neighborhood of BRAF and MAP2K1
biopax_processor = biopax.process_pc_pathsfromto(['BRAF', 'RAF1'], ['MAP2K1', 'MAP2K2'])
```
At this point, `biopax_processor.statements` contains a list of INDRA 
Statements extracted from the paths-from-to query.

## Funding

The development of INDRA has been funded from the following sources:

| Program                                          | Grant number         |
|--------------------------------------------------|----------------------|
| DARPA Big Mechanism                              | W911NF-14-1-0397     |
| DARPA World Modelers                             | W911NF-18-1-0014     |
| DARPA Communicating with Computers               | W911NF-15-1-0544     |
| DARPA Automated Scientific Discovery Fraemwork   | W911NF018-1-0124     |
| DARPA Automating Scientific Knowledge Extraction | HR00111990009        |
