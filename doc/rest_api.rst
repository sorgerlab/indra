REST API
========

Many functionalities of INDRA can be used via a REST API. This enables
making use of INDRA's knowledge sources and assembly capabilities in a
RESTful, platform independent fashion. The REST service is available as a
public web service at http://api.indra.bio:8000 and can also be run locally.

Local installation and use
--------------------------
Running the REST service requires the `flask`, `flask_restx` and `flask_cors`
packages to be installed in addition to all the other requirements of INDRA.
The REST service can be launched by running `api.py` in the `rest_api` folder
within `indra`.

As an alternative, the REST service can be run via the INDRA Docker without
the need for installing any dependencies as follows:

.. code-block:: sh

    docker pull labsyspharm/indra
    docker run -id -p 8080:8080 --entrypoint python labsyspharm/indra /sw/indra/rest_api/api.py

Documentation
-------------
The specific end-points and input/output parameters offered by the REST API
are documented at http://api.indra.bio:8000 or the local address on which
the API is running.
