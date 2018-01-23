from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import pickle
import sys
import logging
import re

logger = logging.getLogger('evidence_subtype_tagger')

reach_rule_filename = 'data/reach_rule_substrings.p'
try:
    reach_rule_regexp = pickle.load( open( reach_rule_filename, "rb" ) )
except:
    logger.error('Cannot open ', reach_rule_filename, ', will not be able' + 
            ' to use priors for individual reach rules')

def determine_reach_subtype(evidence):
    """Returns the subtype of the reach rule. Looks at a list of regular
    expressiosn corresponding to reach rule types, and returns the longest regexp that matches, or None if none of them match.
    
    Parameters
    ----------
    evidence: indra.statements.Evidence
        A reach evidence object to subtype
    """

    assert(evidence.source_api == 'reach')
    event_name = evidence.annotations['found_by']

    best_match_length = None
    best_match = None
    for ss in reach_rule_regexp:
        if re.search(ss, event_name):
            if best_match is None or len(ss) > best_match_length:
                best_match = ss
                best_match_length = len(ss)

    return best_match


def tag_evidence_subtype(evidence):
    """Returns the type and subtype of an evidence object as a string, typically the extraction
    rule or database from which the statement was generated.

    For biopax, this is just the database name.

    Parameters
    ----------
    statement: indra.statements.Evidence
        The statement which we wish to subtype

    Returns
    -------
    types: tuple
        A tuple with (type, subtype), both strings
        Returns (type, None) if the type of statement is not yet handled in
        this function.
    """

    source_api = evidence.source_api
    annotations = evidence.annotations

    if source_api == 'biopax':
        subtype = annotations['source_sub_id']
    elif source_api == 'reach':
        subtype = determine_reach_subtype(evidence)
    elif source_api == 'geneways':
        subtype = annotations['actiontype']
    else:
        subtype = None

    return (source_api, subtype)



