The HTML Curation Interface
===========================
Accessing the Interface
-----------------------
You will usually access this interface from any INDRA application that
exposes statements to you. However if you just want to try out the interface
or don't want to take the detour through any of the applications, you can
follow the format below to access the interface directly in your browser from
the rest API::

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

1) A valid **API key** (at the top of the page, see image)
2) A **curator ID**, such as name or email (at the top of the page, see image)
3) A **selection in the dropdown menu** (by the curated statement)

If you selected *other...* in the dropdown menu, you must *also* describe the
error in the textbox.

+-----------------------------------------+
| .. figure:: images/apikey_curatorID.png |
|   :align: center                        |
|   :figwidth: 75 %                       |
+-----------------------------------------+

When you have entered the necessary information, click the **Submit button** by
the statement that you curated:

+------------------------------------------------+
| .. figure:: images/curation_submit_circled.png |
|   :align: center                               |
|   :figwidth: 75 %                              |
+------------------------------------------------+

A status message will appear once a the server has processed the submission,
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

Curation Best Practices
-----------------------
- Please be consistent in which curator ID you are using. Keeping track of who
  curated what really helps us to faster track down issues with readers and
  the processes that generate statements.
- If you spot multiple levels of errors in a statement - evidence text pair,
  use the most relevant error type in the dropdown menu. E.g. if you see both
  a grounding error and a polarity error, you should pick the grounding
  error since a statement with a grounding error generally would not exist
  if the grounding was correct.
- If you still feel like multiple errors are appropriate for the curation,
  select a new next error from the dropdown menu and make a new submission.
