"""This module implements an API and processor for annotations coming from
hypothes.is. Annotations for a given group are obtained and processed either
into INDRA Statements or into entity grounding annotations.

Two configurable values (either in the INDRA config file or as an environmental
variable) are used. HYPOTHESIS_API_KEY is an API key used to access the
hypothes.is API. HYPOTHESIS_GROUP is an optional configuration used to
select a specific group of annotations on hypothes.is by default.
"""
from .api import *
