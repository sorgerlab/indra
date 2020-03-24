"""This module implements an API and processor for annotations coming from
hypothes.is. Annotations for a given group are obtained and processed either
into INDRA Statements or into entity grounding annotations.

Two configurable values (either in the INDRA config file or as an environmental
variable) are used. HYPOTHESIS_API_KEY is an API key used to access the
hypothes.is API. HYPOTHESIS_GROUP is an optional configuration used to
select a specific group of annotations on hypothes.is by default.

Curation tutorial
------------------
Go to https://web.hypothes.is/ and create an account, and then create
a group in which annotations will be collected. Under Settings, click on
Developer to find the API key. Set his API key in the INDRA config file
under HYPOTHESIS_API_KEY. Optionally, set the group's ID
as HYPOTHESIS_GROUP in the INDRA config file. (Note that both these
values can also be set as environmental variables.) Next, install the
hypothes.is browser plug-in and log in.

Curating Statements
~~~~~~~~~~~~~~~~~~~
To curate text from a website with the intention of creating one or more INDRA
Statements, select some text and create a new annotation using the
hypothes.is browser plug-in. The content of the annotation consists of
one or more lines. The first line should contain one or more English sentences
describing the mechanism(s) that will be represented as an INDRA Statement
(e.g., AMPK activates STAT3) based on the selected text. Each subsequent
line of the annotation is assumed to be a context annotation. These lines
are of the form "<context type>: <context text>" where <context type> can be one
of: Cell type, Cell line, Disease, Organ, Location, Species, and <context text>
is the text describing the context, e.g., lysosome, liver, prostate cancer, etc.

The annotation should also be tagged with `indra` (though by default, if no
tags are given, the processor assumes that the given annotation is an
INDRA Statement annotation).

Curating grounding
~~~~~~~~~~~~~~~~~~
Generally, grounding annotations are only needed if INDRA's current resources
(reading systems, grounding mapping, Gilda, etc.) don't contain a given
synonym for an entity of interest.

With the hypothes.is browser plug-in, select some text on a website that
contains lexical information about an entity or concept of interest.
The conctent of the new annotation can contain one or more lines with identical
syntax as follows:
[text to ground] -> <db_name1>:<db_id1>|<db_name2>:<db_id2>|...
In each case, db_name is a grounding database name space such as HGNC or CHEBI,
and db_id is a value within that namespace such as 1097 or CHEBI:63637.
Example: [AMPK] -> FPLX:AMPK.

The annotation needs to be tagged with `gilda` for the processor to know that
it needs to be interpreted as a grounding annotation.
"""
from .api import *
