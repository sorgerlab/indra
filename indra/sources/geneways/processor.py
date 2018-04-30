"""
This module provides an input processor for information extracted using the
Geneways software suite, converting extraction data in Geneways format into
INDRA statements.

See publication:
Rzhetsky, Andrey, Ivan Iossifov, Tomohiro Koike, Michael Krauthammer, Pauline
Kra, Mitzi Morris, Hong Yu et al. "GeneWays: a system for extracting,
analyzing, visualizing, and integrating molecular pathway data."
Journal of biomedical informatics 37, no. 1 (2004): 43-53.
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from indra.statements import Evidence, Agent
import indra.databases.hgnc_client as hgc
from indra.literature import *
from indra.statements import Complex, Phosphorylation
from indra.sources.geneways.geneways_action_parser import GenewaysActionParser
try:
    from indra.sources.geneways.find_full_text_sentence import FullTextMention
    get_ft_mention = True
except ImportError:
    logger.error('Install the nltk and stemming packages to extract full '
                 'text evidence for Geneways mentions.')
    get_ft_mention = False


logger = logging.getLogger('geneways')


# This will take in an action and action mention and create a single statement
class GenewaysProcessor(object):
    """The GenewaysProcessors converts extracted Geneways action mentions into
    INDRA statements.

    Parameters
    ----------
    search_path : list[str]
        A list of directories in which to search for Geneways data

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA statements converted from Geneways action
        mentions, populated by calling the constructor
    """
    def __init__(self, search_path, get_evidence=True):
        if get_evidence and get_ft_mention:
            self.get_ft_mention = True
        else:
            self.get_ft_mention = False

        # Parse Geneways data. Will give an error if it can't find
        # the Geneways data
        logger.info('Loading Geneways extractions')
        parser = GenewaysActionParser(search_path)
        logger.info('Geneways extractions loaded')
        actions = parser.actions

        # Make a list of statements from the actions
        self.statements = []
        for action in actions:
            for mention in action.action_mentions:
                if mention.negative != '1':
                    new_statement = self.make_statement(action, mention)
                    if new_statement is not None:
                        self.statements.append(new_statement)

    def make_statement(self, action, mention):
        """Makes an INDRA statement from a Geneways action and action mention.

        Parameters
        ----------
        action : GenewaysAction
            The mechanism that the Geneways mention maps to. Note that
            several text mentions can correspond to the same action if they are
            referring to the same relationship - there may be multiple
            Geneways action mentions corresponding to each action.
        mention : GenewaysActionMention
            The Geneways action mention object corresponding to a single
            mention of a mechanism in a specific text. We make a new INDRA
            statement corresponding to each action mention.

        Returns
        -------
        statement : indra.statements.Statement
            An INDRA statement corresponding to the provided Geneways action
            mention, or None if the action mention's type does not map onto
            any INDRA statement type in geneways_action_type_mapper.
        """
        (statement_generator, is_direct) = \
            geneways_action_to_indra_statement_type(mention.actiontype,
                                                    action.plo)

        if statement_generator is None:
            # Geneways statement does not map onto an indra statement
            return None

        # Try to find the full-text sentence
        # Unfortunately, the sentence numbers in the Geneways dataset
        # don't correspond to an obvious sentence segmentation.
        # This code looks for sentences with the subject, object, and verb
        # listed by the Geneways action mention table and only includes
        # it in the evidence if there is exactly one such sentence
        if self.get_ft_mention:
            try:
                content, content_type = get_full_text(mention.pmid, 'pmid')
                if content is not None:
                    ftm = FullTextMention(mention, content)
                    sentences = ftm.find_matching_sentences()
                    if len(sentences) == 1:
                        text = sentences[0]
                    else:
                        text = None
                else:
                    text = None
            except:
                logger.warning('Could not fetch full text for PMID ' +
                               mention.pmid)
        else:
            text = None

        # Make an evidence object
        epistemics = dict()
        epistemics['direct'] = is_direct
        annotations = mention.make_annotation()
        annotations['plo'] = action.plo  # plo only in action table
        evidence = Evidence(source_api='geneways',
                            source_id=mention.actionmentionid,
                            pmid=mention.pmid, text=text,
                            epistemics=epistemics,
                            annotations=annotations)

        # Construct the grounded and name standardized agents
        # Note that this involves grounding the agent by
        # converting the Entrez ID listed in the Geneways data with
        # HGNC and UniProt
        upstream_agent = get_agent(mention.upstream, action.up)
        downstream_agent = get_agent(mention.downstream, action.dn)

        # Make the statement
        return statement_generator(upstream_agent, downstream_agent, evidence)


def get_agent(raw_name, entrez_id):
    db_refs = {'TEXT': raw_name}
    logger.debug('Looking up grounding data for Entrez #%s' % entrez_id)
    hgnc_id = hgc.get_hgnc_from_entrez(entrez_id)
    if hgnc_id is not None:
        db_refs['UP'] = hgc.get_uniprot_id(hgnc_id)
        name = hgc.get_hgnc_name(hgnc_id)
    else:
        name = raw_name
    agent = Agent(name, db_refs=db_refs)
    return agent


def geneways_action_to_indra_statement_type(actiontype, plo):
    """Return INDRA Statement corresponding to Geneways action type.

    Parameters
    ----------
    actiontype : str
        The verb extracted by the Geneways processor
    plo : str
        A one character string designating whether Geneways classifies
        this verb as a physical, logical, or other interaction

    Returns
    -------
    statement_generator :
        If there is no mapping to INDRA statements from this action type
        the return value is None.
        If there is such a mapping, statement_generator is an anonymous
        function that takes in the subject agent, object agent, and evidence,
        in that order, and returns an INDRA statement object.
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
