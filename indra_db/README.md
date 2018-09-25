[![Build Status](https://travis-ci.org/indralab/indra_db.svg?branch=travis_ci)](https://travis-ci.org/indralab/indradb) [![Documentation Status](https://readthedocs.org/projects/indra_db/badge/?version=latest)](https://indradb.readthedocs.io/en/latest/?badge=latest)

# INDRA DB

<img align="left" src="/doc/indra_db_logo.png?raw=True" width="300" height="224" />

The INDRA Database a framwork for creating, maintaining, and accessing a
database (specifically implemented for AWS RDS Postrgres 9+) of content,
readings, and statements. Used as a backend to INDRA, the INDRA Database
provides a systematic way of scaling the knowledge at your fingertips.

### Knowledge sources

The INDRA Database currently integrates the following natural language
processing systems at scale:
- [REACH](https://github.com/clulab/reach) - for biology
- [Sparser](https://github.com/ddmcdonald/sparser) - for biology

with content drawn from:
- [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/) - ~28 million abstracts
- [PubMed Central](/www.ncbi.nlm.nih.gov/pmc/) - ~5 million fulltext
- [Elsevier](https://www.elsevier.com/) - fulltext (requires special access)

We also collect information from these databases:
- [Pathway Commons database](http://pathwaycommons.org/) or any source
    using the [BioPAX](http://www.biopax.org/) format
- [BEL Large Corpus](https://github.com/OpenBEL/) or any source using the
    [BEL](https://github.com/OpenBEL/) format
- [SIGNOR](https://signor.uniroma2.it/)

These databases are retrieved using the tools in `indra.sources`. The statements
extracted from all of these sources are stored and updated in the database.

### Knowledge Assembly

The INDRA Database uses the powerful internal assembly tools available in INDRA
but implemented for large-scale incremental assembly. The resulting corpus of
cleaned and de-duplicated statements, each with fully maintained provenance, is
the primary product of the database.

For more details on the internal assembly process of INDRA, see the
[INDRA documentation](http://indra.readthedocs.io/en/latest/modules/preassembler).

### Access

The content in the database can be accessed by those that created it using the
`indra_db.client` submodule. This repo also implements a REST API which can be
used by those without direct acccess to the database. For access to our REST
API, please contact the authors.

## Installation

The INDRA database only works for Python 3 (tested in 3.5 and 3.6).

First, [install INDRA](http://indra.readthedocs.io/en/latest/installation.html),
then simply clone this repo, and make sure that it is visible in your
`PYTHONPATH`.
