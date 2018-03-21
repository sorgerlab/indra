"""
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

For eidos to provide grounding information to be included in INDRA Statements,
the eidos configuration needs to be adjusted.
First, download vectors.txt from
https://drive.google.com/open?id=1tffQuLB5XtKcq9wlo0n-tsnPJYQF18oS and
put it in a folder called src/main/resources/org/clulab/wm/eidos/w2v within
your eidos folder.
Next, set the property "useW2V" to true in src/main/resources/eidos.conf.
Finally, rerun sbt assembly.
"""

from .eidos_api import process_text, \
    process_json_str, process_json_file, process_json, \
    process_json_ld_str, process_json_ld_file, process_json_ld
