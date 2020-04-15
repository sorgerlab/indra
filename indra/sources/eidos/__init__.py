"""
Eidos is an open-domain machine reading system which uses a cascade of grammars
to extract causal events from free text. It is ideal for modeling applications
that are not specific to a given domain like molecular biology.

To cover a wide range of use cases and scenarios, there are currently 5
different ways in which INDRA can use Eidos.

In all cases for Eidos to provide grounding information to be included in
INDRA Statements, it needs to be configured explicitly to do so.
Please follow instructions at https://github.com/clulab/eidos#configuring to
download and configure Eidos grounding resources.


1. INDRA communicating with a separately running Eidos webapp (:py:mod:`indra.sources.eidos.client`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup and usage: Clone and run the Eidos web server.

.. code-block:: bash

    git clone https://github.com/clulab/eidos.git
    cd eidos
    sbt webapp/run

Then read text by specifying the webserver parameter when using
`indra.sources.eidos.process_text`.

.. code-block:: python

   from indra.sources import eidos
   ep = eidos.process_text('rainfall causes floods',
                            webservice='http://localhost:9000')

Advantages:

* Does not require setting up the pyjnius Python-Java bridge
* Does not require assembling an Eidos JAR file

Disadvantages:

* Not all Eidos functionalities are immediately exposed through its webapp.

2. INDRA using an Eidos JAR directly through a Python-Java bridge (:py:mod:`indra.sources.eidos.reader`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup and usage:

First, the Eidos system and its dependencies need to be packaged as a fat JAR:

.. code-block:: bash

    git clone https://github.com/clulab/eidos.git
    cd eidos
    sbt assembly

This creates a JAR file in eidos/target/scala[version]/eidos-[version].jar.
Set the absolute path to this file on the EIDOSPATH environmental variable
and then append EIDOSPATH to the CLASSPATH environmental variable (entries
are separated by colons).

The `pyjnius` package needs to be set up and be operational. For more details,
see :ref:`pyjniussetup` setup instructions in the documentation.

Then, reading can be done simply using the `indra.sources.eidos.process_text`
function.

.. code-block:: python

   from indra.sources import eidos
   ep = eidos.process_text('rainfall causes floods')

Advantages:

* Doesn't require running a separate process for Eidos and INDRA
* Having a single Eidos JAR file makes this solution portable

Disadvantages:

* Requires configuring pyjnius which is often difficult
* Requires building a large Eidos JAR file which can be time consuming
* The EidosReader instance needs to be instantiated every time a new INDRA
  session is started which is time consuming.

3. INDRA using a Flask sever wrapping an Eidos JAR in a separate process (:py:mod:`indra.sources.eidos.server`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup and usage: Requires building an Eidos JAR and setting up pyjnius --
see above.

First, run the server using


.. code-block:: bash

    python -m indra.sources.eidos.server

Then point to the running server with the webservice parameter when calling
`indra.sources.eidos.process_text`.

.. code-block:: python

   from indra.sources import eidos
   ep = eidos.process_text('rainfall causes floods',
                            webservice='http://localhost:6666')

Advantages:

* EidosReader is instantiated by the Flask server in a separate process,
  therefore it isn't reloaded each time a new INDRA session is started
* Having a single Eidos JAR file makes this solution portable

Disadvantages:

* Currently does not offer any additional functionality compared to running
  the Eidos webapp directly
* Requires configuring pyjnius which is often difficult
* Requires building a large Eidos JAR file which can be time consuming

4. INDRA calling the Eidos CLI using java through the command line (:py:mod:`indra.sources.eidos.cli`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup and usage: Requires building an Eidos JAR and setting EIDOSPATH but does
not require setting up pyjnius -- see above. To use, call any of the
functions exposed in :py:mod:`indra.sources.eidos.cli`.

Advantages:

* Provides a Python-interface for running Eidos on "large scale" jobs, e.g.,
  a large number of input files.
* Does not require setting up pyjnius since it uses Eidos via the command line.
* Provides a way to use any available entrypoint of Eidos.

Disadvantages:

* Requires building an Eidos JAR which can be time consuming.

5. Use Eidos separately to produce output files and then process those with INDRA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this usage mode Eidos is not directly invoked by INDRA. Rather, Eidos
is set up and run idenpendently of INDRA to produce JSON-LD output files
for a set of text content.
One can then use :py:mod:`indra.sources.eidos.api.process_json_file`
in INDRA to process the JSON-LD output files.
"""

from .api import process_text, process_json_str, process_json_file, \
    process_json, reground_texts, initialize_reader, process_json_bio, \
    process_text_bio
