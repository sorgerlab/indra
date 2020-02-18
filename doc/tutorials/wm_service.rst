Summary
-------
These instructions allow you to have reading and assembly services running
locally during the evaluation. With this setup, CauseMos can interact with
the following endpoints:

- http://localhost:8000/reader/process_text to read text via
  eidos/sofia/cwms and get back INDRA Statements
- http://localhost:8001/action to do the following actions: submit_curation,
  reset_curation, update_beliefs, add_ontology_entry, reset_ontology,
  update_groundings.

The instructions below run each Docker container with the `-d` option which
will run containers in the background. You can list running containers with
`docker ps` and stop a container with `docker stop <container id>`.

For interactive text reading, it makes sense to do an initial test reading
before the demo so that Eidos loads the necessary resources. Subsequent
reading calls will be much faster.

Setting up the Eidos service
----------------------------
Clone the Eidos repo and cd to Docker folder

.. code-block:: sh

    git clone https://github.com/clulab/eidos.git
    cd eidos/Docker

Build the Eidos docker image
.. code-block:: sh

    docker build -f Dockerfile . -t eidos-webservice

Run the Eidos web service and expose it on port 9000

.. code-block:: sh

    docker run -id -p 9000:9000 eidos-webservice


Building an Eidos JAR
---------------------
In addition, for the live curation service (in particular, to support
real-time regrounding), a JAR-packaged version of Eidos needs to be
available. In the top-level Eidos folder, run

.. code-block:: sh

    sbt assembly

to obtain a file called something like `target/scala-2.12/eidos-assembly-0.2
.3-SNAPSHOT.jar`. We will refer to this file in later sections as the Eidos
JAR. Note the following details:

- Eidos needs to be configured to use grounding, and one of two large
  word2vec files needs to be made available to Eidos. The commands to do
  this are:

.. code-block:: sh

    wget https://s3.amazonaws.com/world-modelers/data/vectors.txt
    mv vectors.txt src/main/resources/org/clulab/wm/eidos/english/w2v/
    sed -i 's/useW2V = false/useW2V = true/' src/main/resources/eidos.conf
    sed -i 's/useW2V = false/useW2V = true/' src/main/resources/reference.conf
    sed -i 's/glove.840B.300d.txt/vectors.txt/' src/main/resources/eidos.conf


- Running :code:`sbt assembly` may fail with an out of memory error, if this
  happens, run it again with an extra flag to allocate more memory, like
  :code:`sbt -J-Xmx8G assembly`.
- Note that it is also possible to obtain the Eidos JAR from a running
  instance of the Eidos service Docker container using :code:`docker cp`
  thereby avoiding having to build it manually.

Setting up the general INDRA service
------------------------------------
Pull the INDRA docker image from DockerHub

.. code-block:: sh

    docker pull labsyspharm/indra

Run the INDRA web service and expose it on port 8000

.. code-block:: sh

    docker run -id -p 8000:8080 --entrypoint python labsyspharm/indra /sw/indra/rest_api/api.py

Setting up the live feedback service
------------------------------------
Assuming you already have the INDRA docker image and the Eidos JAR, run the
INDRA live feedback service with the following parameters:

- :code:`<folder with corpus>` needs to be a folder in which you have the
  corpus file as well as the Eidos JAR. This folder will be mounted into the
  Docker container on the :code:`/sw/mounted` path allowing the container to
  access both the corpus file, and the Eidos JAR.
- :code:`<corpus file>` needs to be the name of the corpus file in the
  :code:`<folder with corpus>` folder.

.. code-block:: sh

    docker run -v <folder with corpus>:/sw/mounted -id \
      -p 8001:8001 -e EIDOSPATH=/sw/mounted/eidos-assembly-0.2.3-SNAPSHOT.jar \
      --entrypoint python labsyspharm/indra \
      /sw/indra/indra/tools/live_curation.py --json /sw/mounted/<corpus file>

Using the live feedback service
-------------------------------
Below each example uses the remote service, you can replace that IP with
localhost to do the same locally.

Submit curations for a set of Statements in a corpus:

.. code-block::

    URL: http://54.84.114.146:8001/submit_curation
    Method: POST with JSON content header
    Input parameters: {"corpus_id": "1", "curations": {"38ce0c14-2c7e-4df8-bd53-3006afeaa193": 0}}
    Output: {}

Update beliefs of a corpus:

.. code-block::

    URL: http://54.84.114.146:8001/update_beliefs
    Method: POST with JSON content header
    Input parameters: {"corpus_id": "1"}
    Output: {"38ce0c14-2c7e-4df8-bd53-3006afeaa193": 0,
     "6f2b2d69-16af-40ea-aa03-9b3a9a1d2ac3": 0.6979166666666666,
     "727adb95-4890-4bbc-a985-fd985c355215": 0.6979166666666666}


Reset all submitted curations so far:

.. code-block::

    URL: http://54.84.114.146:8001/reset_curation
    Method: POST with JSON content header
    Input parameters: {}
    Output: {}

Add a new ontology entry:

.. code-block::

    URL: http://54.84.114.146:8001/add_ontology_entry
    Method: POST with JSON content header
    Input parameters: {"entry": "UN/animals/dog", "examples": ["dog", "canine", "puppy"]}
    Output: {}

Reset all customizations to the ontology so far:

.. code-block::

    URL: http://54.84.114.146:8001/reset_ontology
    Method: POST with JSON content header
    Input parameters: {}
    Output: {}

Update groundings and re-assemble corpus based on current ontology:

.. code-block::

    URL: http://54.84.114.146:8001/update_groundings
    Method: POST with JSON content header
    Input parameters: {"corpus_id": "1"}
    Output: [{"type": "Influence", ...}] (INDRA Statements JSON)
