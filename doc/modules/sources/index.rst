Processors for knowledge input (:py:mod:`indra.sources`)
========================================================

INDRA interfaces with and draws knowledge from many sources including
reading systems (some that extract biological mechanisms, and some that extract
general causal interactions from text) and also from structured databases,
which are typically human-curated or derived from experimental data.

Reading Systems
---------------
.. toctree::
   :maxdepth: 3

   reach/index
   trips/index
   sparser/index
   medscan/index
   tees/index
   isi/index
   geneways/index
   rlimsp/index
   eidos/index

Molecular Pathway Databases
---------------------------
.. toctree::
   :maxdepth: 3

   bel/index
   biopax/index
   signor/index
   biogrid/index
   hprd/index
   trrust
   phosphoelm/index
   virhostnet/index
   omnipath/index

Chemical Information Databases
------------------------------
.. toctree::
   :maxdepth: 3

   ctd/index
   drugbank/index
   dgi/index
   tas/index
   crog/index

Custom Knowledge Bases
----------------------
.. toctree::
   :maxdepth: 3

   ndex_cx/index
   indra_db_rest/index
   hypothesis/index
   biofactoid/index
   minerva/index

Utilities
---------

.. automodule:: indra.sources.utils
    :members:
