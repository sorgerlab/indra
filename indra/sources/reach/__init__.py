"""
REACH is a biology-oriented machine reading system which uses a cascade of
grammars to extract biological mechanisms from free text.

To cover a wide range of use cases and scenarios, there are currently 4
different ways in which INDRA can use REACH.


1. INDRA communicating with a remote REACH API Server (:py:mod:`indra.sources.reach.api`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup and usage: Does not require any additional setup after installing INDRA.

Read text using the default values for `offline` and `url` parameters.

.. code-block:: python

   from indra.sources import reach
   rp = reach.process_text('MEK binds ERK')

It is also possible to read NXML (string or file) and process text of a paper
given its PMC ID or PubMed ID.

Advantages:

* Does not require setting up the pyjnius Python-Java bridge
* Does not require assembling a REACH JAR file

Disadvantages:

* No control on which REACH version is used to run the service
* Difficulties processing NXML-formatted text have been observed in the past

2. INDRA communicating with a locally running REACH API Server (:py:mod:`indra.sources.reach.api`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup and usage: Clone and run REACH web server.

.. code-block:: bash

    git clone https://github.com/clulab/reach.git
    cd reach
    sbt 'run-main org.clulab.reach.export.server.ApiServer'

Then read text by specifying the url parameter when using
`indra.sources.reach.process_text`.

.. code-block:: python

   from indra.sources import reach
   rp = reach.process_text('MEK binds ERK',
                           url='http://localhost:8080/api/text')

It is also possible to read NXML (string or file) and process text of a paper
given its PMC ID or PubMed ID.

Advantages:

* Control the REACH version used to run the service
* Does not require setting up the pyjnius Python-Java bridge
* Does not require assembling a REACH JAR file

Disadvantages:

* First request might be time-consuming as REACH is loading additional
resources.

3. INDRA using a REACH JAR directly through a Python-Java bridge (:py:mod:`indra.sources.reach.reader`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup and usage:

First, the REACH system and its dependencies need to be packaged as a fat JAR:

.. code-block:: bash

    git clone https://github.com/clulab/reach.git
    cd reach
    sbt assembly

This creates a JAR file in reach/target/scala[version]/reach-[version].jar.
Set the absolute path to this file on the REACHPATH environmental variable
and then append REACHPATH to the CLASSPATH environmental variable (entries
are separated by colons).

The `pyjnius` package needs to be set up and be operational. For more details,
see :ref:`pyjniussetup` setup instructions in the documentation.

Then, reading can be done simply using the `indra.sources.reach.process_text`
function with offline option.

.. code-block:: python

   from indra.sources import reach
   rp = reach.process_text('MEK binds ERK', offline=True)


Advantages:

* Doesn't require running a separate process for REACH and INDRA
* Having a single REACH JAR file makes this solution portable

Disadvantages:

* Requires configuring pyjnius which is often difficult
* Requires building a large REACH JAR file which can be time consuming
* The ReachReader instance needs to be instantiated every time a new INDRA
  session is started which is time consuming.

4. Use REACH separately to produce output files and then process those with INDRA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this usage mode REACH is not directly invoked by INDRA. Rather, REACH
is set up and run idenpendently of INDRA to produce JSON-LD output files
for a set of text content.
One can then use :py:mod:`indra.sources.reach.api.process_json_file`
in INDRA to process the JSON-LD output files.
"""


from .api import (process_pmc,
                  process_pubmed_abstract,
                  process_text,
                  process_nxml_str,
                  process_nxml_file,
                  process_json_str,
                  process_json_file)
