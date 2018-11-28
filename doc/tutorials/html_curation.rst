The Statement curation interface
================================

You will usually access this interface from an INDRA application that
exposes statements to you. However if you just want to try out the interface
or don't want to take the detour through any of the applications, you can
follow the format below to access the interface directly in your browser from
the INDRA-DB REST API::

    http://api.host/statements/from_agents?subject=SUBJ&object=OBJ&api_key=12345&format=html

where *api.host* should be replaced with the address to the REST API service
(see the `documentation
<https://github.com/indralab/indra_db/blob/master/rest_api/README.md>`_).
Entering the whole address in your browser will query for statements where
*SUBJ* is the subject and *OBJ* is the object of the statements.

For more details about what options are available when doing curation, please
refer to the `curation section
<https://github.com/indralab/indra_db/blob/master/rest_api/README.md#curation>`_
of the documentation.

Curating a Statement
--------------------
Let's assume you want to check any statements were ROS1 is an agent for
errors. Let's also limit the number of statements to 100 and the number of
evidences per statements to 10. This will speed up the query and page loading.
The appropriate address to enter in your browser would then be::

    http://api.host/statements/from_agents?agent=ROS1&format=html&ev_limit=10&max_stmts=100

To start curating a statement, **click the pen icon (circled)** on the far left
side of the statement. This will produce a row below the statement with a
dropdown menu, a text box and a submit button:

+-----------------------------------------------------+
| .. figure:: images/curation_row_created_circled.png |
|   :align: center                                    |
|   :figwidth: 75 %                                   |
+-----------------------------------------------------+

The **dropdown menu** contains common errors and also the possibility to mark
the statement as 'correct'. If none of the types fit, select the *other...*
option, and describe the error with one or a few words in the provided
textbox. Note that if you pick *other...*, describing the error is mandatory.
In our example, we see that *reactive oxygen species* is incorrectly grounded
to *ROS*, so we pick *grounding* from the dropdown menu:

+------------------------------------------------------+
| .. figure:: images/curation_select_error_circled.png |
|    :align: center                                    |
|    :figwidth: 75 %                                   |
+------------------------------------------------------+

In the textbox, you can add a short optional description to clarify why you
marked this piece of evidence with the error type you chose. When you are
done, you are ready to submit your curation.

Submitting a Curation
---------------------
To **submit a curation**, there are three minimum requirements:

1) A valid **API key** (at the top of the page, see image below)
2) A **curator ID**, such as name or email (at the top of the page, see image
   below)
3) A **selection in the dropdown menu** (by the curated statement)

If you selected *other...* in the dropdown menu, you must *also* describe the
error in the textbox.

+------------------------------------------+
| .. figure:: images/apikey_curatorID.png  |
|   :align: center                         |
|   :figwidth: 75 %                        |
|                                          |
|   *Provide a valid API key and a user ID |
|   in order to submit a curation*         |
+------------------------------------------+

When you have entered the necessary information, click the **Submit button** by
the statement that you curated:

+------------------------------------------------+
| .. figure:: images/curation_submit_circled.png |
|   :align: center                               |
|   :figwidth: 75 %                              |
+------------------------------------------------+

A status message will appear once the server has processed the submission,
indicating if the submission was successful or which problem arose if not.
The pen icon will also change color based in the returned status. **Green**
indicates a successful submission:

+--------------------------------------------------------------+
| .. figure:: images/curation_submitted_successfully.png       |
|   :align: center                                             |
|   :figwidth: 75 %                                            |
|                                                              |
|   *A green icon indicates a successfully submitted curation* |
+--------------------------------------------------------------+

while a **red** indicates something went wrong with the submission:

+--------------------------------------------------------------------------+
| .. figure:: images/bad_submission.png                                    |
|   :align: center                                                         |
|   :figwidth: 75 %                                                        |
|                                                                          |
|   *A red icon indicates that something went wrong during the submission* |
+--------------------------------------------------------------------------+

Curation Guidelines
-------------------
Basic principles
~~~~~~~~~~~~~~~~
The main question to ask when deciding whether a given Statement is correct
with respect to a given piece of evidence is::

    Is there support in the evidence sentence for the Statement?

If the answer is **Yes**, then the given sentence
is a valid piece of evidence for the Statement. In fact, you can assert this
correctness by choosing the "Correct" option from the curation drop-down list.
Curations that assert correctness are just as valuable as curations of
incorrectness so the use of this option is encouraged.

Assuming the answer to the above question is **No**, one needs to determine
what the error can be attributed to. The following section describes the
specific error types that can be flagged.

Types of errors to curate
~~~~~~~~~~~~~~~~~~~~~~~~~
There are currently the following options to choose from when curating
incorrect Statement-sentence relationships:

* **Entity Boundaries**: this is applicable if the bounderies of one of the named
  entities was incorrectly recognized. Example: "gap" is highlighted as an
  entity, when in fact, the entity mentioned in the sentence was
  "gap junction". These errors in entity boundaries almost always result in
  incorrect grounding, since the wrong string is attempted to be grounded.
  Therefore this error "subsumes" grounding errors.
  Note: to help correct entity boundaries, add the following to the
  Optional description text box: [gap junction], i.e. the desired entity
  name inside square brackets.
* **Grounding**: this is applicable if a named entity is assigned an incorrect
  database identifier. Example::

    Assume that in a sentence, "ER" is mentioned referring to endoplasmic
    reticulum, but in a Statement extracted from the sentence, it is
    grounded to the ESR1 (estrogen receptor alpha) gene.

  Note: to help correct grounding, add the following to the Optional
  description text box::

    [ER] -> MESH:D004721

  where [ER] is the entity string,
  MESH is the namespace of a database/ontology, and D004721 is the unique ID
  corresponding to endoplasmic reticulum in MESH.
  A list of commonly used namespaces in INDRA are given in:
  https://indra.readthedocs.io/en/latest/modules/statements.html.
  Note that you can also add multiple groundings separated by "|", e.g.
  HGNC:11998|UP:P04637.

* **Polarity**: this is applicable if an essentially correct Statement was
  extracted but the Statement has the wrong polarity, e.g. Activation
  instead of Inhibition, of Phosphorylation instead of Dephosphorylation.
  Example::

    Sentence: "NDRG2 overexpression specifically inhibits SOCS1 phosphorylation"
    Statement: Phosphorylation(NDRG2(), SOCS1())

  has incorrect polarity. It should be Dephosphorylation instead of
  Phosphorylation.

* **No Relation**: this is applicable if the sentence does not imply a
  relationship between the agents appearing in the Statement. Example::

    Sentence: "Furthermore, triptolide mediated inhibition of NF-kappaB
        activation, Stat3 phosphorylation and increase of SOCS1 expression in
        DC may be involved in the inhibitory effect of triptolide."
    Statement: Phosphorylation(STAT3(), SOCS1())

  can be flagged as No Relation.

* **Wrong Relation Type**: this is applicable if the sentence implies a
  relationship between agents appearing in the Statement but the type of
  Statement is inconsistent with the sentence. Example::

    Sentence: "We report the interaction between tacrolimus and chloramphenicol
        in a renal transplant recipient."
    Statement: Complex(tacrolimus(), chloramphenicol())

  can be flagged as Wrong Relation Type since the sentence implies a drug
  interaction that does not involve complex formation.

* **Activity vs. Amount**: this is applicable when the sentence implies a
  regulation of amount but the corresponding Statement implies regulation
  of activity or vice versa. Example::

    Sentence: "NFAT upregulates IL2"
    Sentence: Activation(NFAT(), IL2())

  Here the sentence implies upregulation of the amount of IL2 but the
  corresponding Statement is of type Activation rather than IncreaseAmount.

* **Negative Result**: this is applicable if the sentence implies the lack of
  or opposite of a relationship. Example::

    Sentence: "These results indicate that CRAF, but not BRAF phosphorylates
        MEK in NRAS mutated cells."
    Statement: Phosphorylation(BRAF(), MEK())

  Here the sentence does not support the Statement due to a negation and
  should therefore be flagged as a Negative Result.

* **Hypothesis**: this is applicable if the sentence describes a hypothesis or
  an experiment rather than a result or mechanism. Example::

    Sentence: "We tested whether EGFR activates ERK."
    Statement: Activation(EGFR(), ERK())

  Here the sentence describes a hypothesis with respect to the Statement, and
  should therefore be flagged as a Hypothesis upon curation (unless of course
  the Statement already has a correct *hypothesis* flag).

* **Agent Conditions**: this is applicable if one of the Agents in the Statement
  is missing relevant conditions that are mentioned in the sentence, or has
  incorrect conditions attached to it. Example::

    Sentence: "Mutant BRAF activates MEK"
    Statement: Activation(BRAF(), MEK())

  can be curated to be missing Agent conditions since the mutation on BRAF is
  not captured.

* **Modification Site**: this is applicable if an amino-acid site is
  missing or incorrect in a modification Statement. Example::

    Sentence: "MAP2K1 phosphorylates MAPK1 at T185."
    Statement: Phosphorylation(MAP2K1(), MAPK1())

  Here the obvious modification site is missing from MAPK1.

* **Other**: this is an option you can choose whenever the problem isn't
  well captured by any of the more specific options. In this case you need
  to add a note to explain what the issue is.

General notes on curation
~~~~~~~~~~~~~~~~~~~~~~~~~

* If you spot multiple levels of errors in a Statement-sentence pair,
  use the most relevant error type in the dropdown menu. E.g. if you see both
  a grounding error and a polarity error, you should pick the grounding
  error since a statement with a grounding error generally would not exist
  if the grounding was correct.
* If you still feel like multiple errors are appropriate for the curation,
  select a new error from the dropdown menu and make a new submission.
* Please be consistent in using your email address as your curator ID.
  Keeping track of who curated what helps us to faster track down
  issues with readers and the assembly processes that generate statements.
