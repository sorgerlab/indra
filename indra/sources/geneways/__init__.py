"""An input processor for information extracted via the Geneways software suite

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
from indra.statements import Evidence, Agent
import indra.databases.hgnc_client as hgc
import logging

logger = logging.getLogger('geneways')


# This will take in an action and action mention and create a single statement
class GenewaysProcessor(object):
    def __init__(self):
        """Parses the Geneways data and makes a list of indra statements"""
        # Parse Geneways data. Will give an error if it can't find
        # the Geneways data
        parser = GenewaysActionParser()
        actions = parser.actions

        # Make a list of statements from the actions
        self.statements = []
        for action in actions:
            new_statement = self.make_statement(action)
            if new_statement is not None:
                self.statements.append(new_statement)

    def make_statement(self, action):
        """Given a single parsed Geneways action, creates an Indra statement.
        Returns None if the Geneways action type is not mapped onto an
        Indra statement type."""
        (statement_generator, is_direct) = \
                geneways_action_to_indra_statement_type(action.actiontype, \
                action.plo)

        if statement_generator is None:
            # Geneways statement does not map onto an indra statement
            return None 
        else:
            # Make an evidence object
            epistemics = dict()
            epistemics['direct'] = is_direct
            pmids = list()
            for mention in action.action_mentions:
                pmids.append(mention.pmid)
            annotations = action.make_annotation()
            evidence = Evidence(source_api='geneways',
                                source_id=action.hiid,
                                pmid=pmids, text=None,
                                epistemics=epistemics,
                                annotations=annotations)

            # Ground the upstream agent
            # Note that we are using the name as it appeared in the text, rather
            # than some standardized name in a database, but grounding it by
            # converting the Entrez ID listed in the Geneways data with
            # HGNC and UniProt
            up_name = action.action_mentions[0].upstream 
            upstream_db = dict()
            logger.debug('Looking up grounding data for Entrez #%s' % action.up)
            upstream_db['HGNC'] = hgc.get_hgnc_from_entrez(action.up)
            upstream_db['UP']   = hgc.get_uniprot_id(upstream_db['HGNC'])
            upstream_db['TEXT']   = up_name
            upstream_agent = Agent(up_name, db_refs=upstream_db)

            # Ground the downstream agent
            down_name = action.action_mentions[0].downstream
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


