World Modelers INDRA service stack
==================================

Setting up the services locally
-------------------------------
These instructions describe setting up and using the INDRA service stack
for World Modelers applications, in particular, as a back-end for the
CauseMos UI.

The instructions below run each Docker container with the :code:`-d` option
which will run containers in the background. You can list running containers
with their ids using :code:`docker ps` and stop a container with
:code:`docker stop <container id>`.

Setting up the Eidos service
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone the Eidos repo and cd to the Docker folder

.. code-block:: sh

    git clone https://github.com/clulab/eidos.git
    cd eidos/Docker

Build the Eidos docker image

.. code-block:: sh

    docker build -f DockerfileRunProd . -t eidos-webservice

Run the Eidos web service and expose it on port 9000

.. code-block:: sh

    docker run -id -p 9000:9000 eidos-webservice


Setting up the general INDRA service
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Pull the INDRA docker image from DockerHub

.. code-block:: sh

    docker pull labsyspharm/indra

Run the INDRA web service and expose it on port 8000

.. code-block:: sh

    docker run -id -p 8000:8080 --entrypoint gunicorn labsyspharm/indra:latest \
        -w 1 -b :8000 -t 600 rest_api.api:app

Note that the :code:`-w 1` parameter specifies one service worker which can
be set to a higher number if needed.

Setting up the INDRA live curation service
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Assuming you already have the INDRA docker image, run the INDRA live
feedback service with the following parameters:

.. code-block:: sh

    docker run -id -p 8001:8001 --env-file docker_variables --entrypoint \
    python labsyspharm/indra /sw/indra/indra/tools/live_curation.py

Here we use the tag :code:`--env-file` to provide a file containing
environment variables to the docker. In this case, we need to provide
:code:`AWS_ACCESS_KEY_ID` and :code:`AWS_SECRET_ACCESS_KEY` to allow the
curation service to access World Modelers corpora on S3.
The file content should look like this:

.. code-block:: sh

    AWS_ACCESS_KEY_ID=<aws_access_key_id>
    AWS_SECRET_ACCESS_KEY=<aws_secret_access_key>

Replace :code:`<aws_access_key_id>` and :code:`<aws_secret_access_key>` with
your aws access and secret keys.

Using the services
------------------
Below, SERVICE_HOST should be replaced by the address of the server on which
the services are running.

Read a given text with a reader and return INDRA Statements (below, <reader>
can be eidos, sofia or cwms). Note that for `eidos` specifically, a
`webservice` parameter should also be passed which points to the address
on which the Eidos web service is running (see above):

.. code-block:: sh

    URL: http://SERVICE_HOST:8000/<reader>/process_text
    Method: POST with JSON content header
    Input parameters: {"text": "rainfall causes floods"}
    Output: {}

Submit curations for a set of Statements in a corpus:

.. code-block:: sh

    URL: http://SERVICE_HOST:8001/submit_curation
    Method: POST with JSON content header
    Input parameters: {"corpus_id": "<corpus-id>", "curations": {"38ce0c14-2c7e-4df8-bd53-3006afeaa193": 0}}
    Output: {}

Save curations for a given corpus on S3:

.. code-block:: sh

    URL: http://SERVICE_HOST:8001/save_curation
    Method: POST with JSON content header
    Input parameters: {"corpus_id": "<corpus-id>"}
    Output: {}

Update beliefs of a corpus:

.. code-block:: sh

    URL: http://SERVICE_HOST:8001/update_beliefs
    Method: POST with JSON content header
    Input parameters: {"corpus_id": "<corpus-id>"}
    Output: {"38ce0c14-2c7e-4df8-bd53-3006afeaa193": 0,
     "6f2b2d69-16af-40ea-aa03-9b3a9a1d2ac3": 0.6979166666666666,
     "727adb95-4890-4bbc-a985-fd985c355215": 0.6979166666666666}

Reset all submitted curations so far:

.. code-block:: sh

    URL: http://SERVICE_HOST:8001/reset_curation
    Method: POST with JSON content header
    Input parameters: {}
    Output: {}

Add a new ontology entry:

.. code-block:: sh

    URL: http://SERVICE_HOST:8001/add_ontology_entry
    Method: POST with JSON content header
    Input parameters: {"entry": "UN/animals/dog", "examples": ["dog", "canine", "puppy"]}
    Output: {}

Reset all customizations to the ontology so far:

.. code-block:: sh

    URL: http://SERVICE_HOST:8001/reset_ontology
    Method: POST with JSON content header
    Input parameters: {}
    Output: {}

Update groundings and re-assemble corpus based on current ontology:

.. code-block:: sh

    URL: http://SERVICE_HOST:8001/update_groundings
    Method: POST with JSON content header
    Input parameters: {"corpus_id": "1"}
    Output: [{"type": "Influence", ...}] (INDRA Statements JSON)
