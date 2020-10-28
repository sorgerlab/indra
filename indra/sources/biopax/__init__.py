"""This module allows processing BioPAX content into INDRA Statements. It uses
the pybiopax package (https://github.com/indralab/pybiopax) to process OWL
files or strings, or to obtain BioPAX content by querying the PathwayCommons
web service. The module has been tested with BioPAX content from PathwayCommons
https://www.pathwaycommons.org/archives/PC2/v12/. BioPAX from other sources may
not adhere to the same conventions and could result in processing issues,
though these can typically be addressed with minor changes in the processor's
logic.
"""
from .api import *
