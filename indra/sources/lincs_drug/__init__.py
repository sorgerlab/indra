"""This module provides and API and processor for the HMS LINCS small molecule
target relationship database. This is a manually curated set of relationships
with the "nominal" target of each drug determined by a human expert. Note that
the determination of the "nominal" target is not always backed up by
experimentally measured affinities. The underlying data is available here:
http://lincs.hms.harvard.edu/db/datasets/20000/results
"""

from .api import *
