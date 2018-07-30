INDRA documentation
===================

INDRA (the Integrated Network and Dynamical Reasoning Assembler) assembles
information about causal mechanisms into a common format that can be used
to build several different kinds of predictive and explanatory models.
INDRA was originally developed for molecular systems biology and is
currently being generalized to other domains.

In molecular biology, sources of mechanistic
information include pathway databases, natural language descriptions of
mechanisms by human curators, and findings extracted from the literature by
text mining.

Mechanistic information from multiple sources is de-duplicated,
standardized and assembled into sets of *Statements* with
associated evidence. Sets of Statements can then be used to assemble both
executable rule-based models (using `PySB`_) and a variety of different types
of network models.

.. _PySB: http://pysb.org

.. toctree::
   :maxdepth: 3

   license.rst
   installation.rst
   getting_started.rst
   modules/index
   tutorials/index
   rest_api.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

