Eidos (:py:mod:`indra.sources.eidos`)
=====================================
Eidos is an open-domain machine reading system which uses a cascade of grammars
to extract causal events from free text. It is ideal for modeling applications
that are not specific to a given domain like molecular biology.

To set up reading with Eidos, the Eidos system and its dependencies need to be
compiled and packaged as a fat JAR:

.. code-block:: bash

    git clone https://github.com/clulab/eidos.git
    cd eidos
    sbt assembly

This creates a JAR file in eidos/target/scala[version]/eidos-[version].jar.
Set the absolute path to this file on the EIDOSPATH environmental variable
and then append EIDOSPATH to the CLASSPATH environmental variable (entries
are separated by colons).

The `pyjnius` package needs to be set up and operational to use Eidos reading
in Python. For more details, see :ref:`pyjniussetup` setup instructions in
the documentation.


Eidos API (:py:mod:`indra.sources.eidos.eidos_api`)
---------------------------------------------------

.. automodule:: indra.sources.eidos.eidos_api
    :members:

Eidos Processor (:py:mod:`indra.sources.eidos.processor`)
---------------------------------------------------------

.. automodule:: indra.sources.eidos.processor
    :members:

Eidos Reader (:py:mod:`indra.sources.eidos.eidos_reader`)
---------------------------------------------------------

.. automodule:: indra.sources.eidos.eidos_reader
    :members:
