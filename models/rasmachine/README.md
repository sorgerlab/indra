Starting a new model
====================

To start a new model, run 

    ./new_machine.sh model_name

This script generates the following folders and files

    model_name
    model_name/log.txt
    model_name/err.txt
    model_name/search_terms.txt
    model_name/belief.txt
    model_name/jsons/
    model_name/jsons/abstract
    model_name/jsons/full

If you want to use Twitter and Gmail bindings, you can put the corresponding
credentials in the model folder.

    model_name/twitter_cred.txt
    model_name/gmail_cred.txt

Extending a model
=================

First, install the machine-specific dependencies::

    pip install -r requirements.txt

To extend a model, run

    ./run_machine.sh model_name

This script sets the JAVA\_HOME, CLASSPATH and PYTHONPATH environmental
variables and then calls

    python rasmachine.py --model model_name [--twitter twitter_cred]
                         [--gmail gmail_cred] [--belief belief_file]

Extending a model involves extracting PMIDs from emails (if Gmail credentials
are given), and searching using INDRA's PubMed client with each line of
search\_terms.txt as a search term.  INDRA's literature client is then used to
find the full text corresponding to each PMID or its abstract when the full
text is not available. The REACH parser is then used to read each new paper.
INDRA uses the REACH output to construct Statements corresponding to
mechanisms.  It then adds them to an incremental model through a process of
assembly involving duplication and overlap resolution and the application of
filters.
