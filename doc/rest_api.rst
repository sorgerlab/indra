REST API
========

Many functionalities of INDRA can be used via a REST API. This enables
making use of INDRA's knowledge sources and assembly capabilities in a
RESTful, platform independent fashion. The REST service is available as a
public web service and can also be run locally (on a single machine or local
network).

Installation
------------
Running the REST service requires the `bottle` package to be installed in
addition to all the other requirements of INDRA.

Launching the REST service (locally)
------------------------------------
The REST service can be launched by running `api.py` in the `rest_api` folder
within `indra`.

Documentation
-------------
The specific end-points and input/output parameters offered by the REST API
are documented in `rest_api/docs/index.html`, which is accessible locally
within the `indra` folder.
