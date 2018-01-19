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

from indra.sources.geneways.geneways_action_parser import GenewaysActionParser
from indra.sources.geneways.geneways_action_type_mapper import \
        geneways_action_to_indra_statement_type
from indra.sources.geneways.find_full_text_sentence import FullTextMention
from indra.statements import Evidence, Agent
import indra.databases.hgnc_client as hgc
import logging
import sys
from indra.literature import *


logger = logging.getLogger('geneways')
logger.setLevel(logging.DEBUG)

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
    def __init__(self, search_path):
        # Parse Geneways data. Will give an error if it can't find
        # the Geneways data
        if sys.version_info[0] < 3:
            logger.warning('This processor is very slow in python 2! ' + 
                    'Python 3 is recommended.')


        logger.debug('Loading Geneways extractions')
        parser = GenewaysActionParser(search_path)
        logger.debug('\tGeneways extractions loaded!')
        actions = parser.actions

        
        # Make a list of statements from the actions
        self.statements = []
        for action in actions:
            for mention in action.action_mentions:
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
            The Geneways action mention object corresponding to a single mention
            of a mechanism in a specific text. We make a new INDRA statement
            corresponding to each action mention.

        Returns
        -------
        statement : indra.statements.Statement
            An INDRA statement corresponding to the provided Geneways action
            mention, or None if the action mention's type does not map onto
            any INDRA statement type in geneways_action_type_mapper.
        """
        (statement_generator, is_direct) = \
                geneways_action_to_indra_statement_type(mention.actiontype, \
                action.plo)

        if statement_generator is None:
            # Geneways statement does not map onto an indra statement
            return None 
        else:
            #Try to find the full-text sentence
            #Unfortunately, the sentence numbers in the Geneways dataset
            #don't correspond to an obvious sentence segmentation.
            #This code looks for sentences with the subject, object, and verb
            #listed by the Geneways action mention table and only includes
            #it in the evidence if there is exactly one such sentence
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
                logger.warning('Could not fetch full text for PMID ' + pmid)


            # Make an evidence object
            epistemics = dict()
            epistemics['direct'] = is_direct
            annotations = mention.make_annotation()
            annotations['plo'] = action.plo #plo only in action table
            evidence = Evidence(source_api='geneways',
                                source_id=mention.actionmentionid,
                                pmid=mention.pmid, text=text,
                                epistemics=epistemics,
                                annotations=annotations)

            # Ground the upstream agent
            # Note that we are using the name as it appeared in the text, rather
            # than some standardized name in a database, but grounding it by
            # converting the Entrez ID listed in the Geneways data with
            # HGNC and UniProt
            up_name = mention.upstream 
            upstream_db = dict()
            logger.debug('Looking up grounding data for Entrez #%s' % action.up)
            upstream_db['HGNC'] = hgc.get_hgnc_from_entrez(action.up)
            upstream_db['UP']   = hgc.get_uniprot_id(upstream_db['HGNC'])
            upstream_db['TEXT']   = up_name
            upstream_agent = Agent(up_name, db_refs=upstream_db)


            # Ground the downstream agent
            down_name = mention.downstream
            downstream_db = dict()
            logger.debug('Looking up grounding data for Entrez #%s' % action.dn)
            downstream_db['HGNC'] = hgc.get_hgnc_from_entrez(action.dn)
            downstream_db['UP']   = hgc.get_uniprot_id(downstream_db['HGNC'])
            downstream_db['TEXT']   = down_name
            downstream_agent = Agent(down_name, db_refs=downstream_db)

            # Make the statement
            return statement_generator(upstream_agent,
                    downstream_agent,
                    evidence)


