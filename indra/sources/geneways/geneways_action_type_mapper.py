from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.statements import Complex, Phosphorylation


def geneways_action_to_indra_statement_type(actiontype, plo):
    """This function provides the mapping from Geneways statements
    to Indra statements. It is used by GenewaysProcessor.
    Inputs:
    * actiontype: the verb extracted by the Geneways processor
    * plo: A one character string designating whether Geneways classifies
    this verb as a physical, logical, or other interaction
    
    Outputs:
    * statement_generator: If there is no mapping to Indra statements
    from this action type, None.
        If there is such a mapping, statement_generator is an anonymous
        function that takes in the subject agent, object agent, and evidence,
        in that order, and returns an Indra statement object.
    """
    actiontype = actiontype.lower()

    statement_generator = None
    is_direct = (plo == 'P')

    if actiontype == 'bind':
        statement_generator = lambda substance1, substance2, evidence: \
                Complex([substance1, substance2], evidence=evidence) 
        is_direct = True
    elif actiontype == 'phosphorylate':
        statement_generator = lambda substance1, substance2, evidence: \
                Phosphorylation(substance1, substance2, evidence=evidence)
        is_direct = True

    return (statement_generator, is_direct)

