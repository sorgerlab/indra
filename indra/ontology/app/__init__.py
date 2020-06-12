"""This module implements IndraOntology functionalities as a web service.
If instantiating an ontology directly is not desirable (for instance because
of memory constraints), this app can be started on a suitable server, and
an instance of the VirtualOntology class can be used to communicate with it
transparently.

To start the server, run

.. code-block:: bash

    python -m indra.ontology.app.app

or use a WSGI application server such as gunicorn (the service uses port 8002
by default, this can be changed using the `--port` argument).

Once the service is
started, one option is to create an instance of
`VirtualOntology(url=<service url>)` and use it as an argument in various
function calls.

Another option is to set the value `INDRA_ONTOLOGY_URL=<service url>` either as
an environmental variable or in the INDRA configuration file. If this value
is set, INDRA will use an appropriate instance of a VirtualOntology
which communicates with the service in place of the BioOntology.
"""
