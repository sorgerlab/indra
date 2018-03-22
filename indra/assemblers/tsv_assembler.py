from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from copy import copy
from indra.databases import get_identifiers_url
from indra.statements import *
from indra.util import write_unicode_csv


logger = logging.getLogger('tsv_assembler')


class TsvAssembler(object):
    """Assembles Statements into a set of tabular files for export or curation.

    Currently designed for use with "raw" Statements, i.e., Statements with a
    single evidence entry. Exports Statements into a single tab-separated file
    with the following columns:

    *INDEX*
        A 1-indexed integer identifying the statement.
    *UUID*
        The UUID of the Statement.
    *TYPE*
        Statement type, given by the name of the class in indra.statements.
    *STR*
        String representation of the Statement. Contains most relevant
        information for curation including any additional statement data
        beyond the Statement type and Agents.
    *AG_A_TEXT*
        For Statements extracted from text, the text in the sentence
        corresponding to the first agent (i.e., the 'TEXT' entry in the
        db_refs dictionary). For all other Statements, the Agent name is
        given. Empty field if the Agent is None.
    *AG_A_LINKS*
        Groundings for the first agent given as a comma-separated list of
        identifiers.org links. Empty if the Agent is None.
    *AG_A_STR*
        String representation of the first agent, including additional
        agent context (e.g. modification, mutation, location, and bound
        conditions). Empty if the Agent is None.
    *AG_B_TEXT, AG_B_LINKS, AG_B_STR*
        As above for the second agent. Note that the Agent may be None (and
        these fields left empty) if the Statement consists only of a single
        Agent (e.g., SelfModification, ActiveForm, or Translocation statement).
    *PMID*
        PMID of the first entry in the evidence list for the Statement.
    *TEXT*
        Evidence text for the Statement.
    *IS_HYP*
        Whether the Statement represents a "hypothesis", as flagged by some
        reading systems and recorded in the `evidence.epistemics['hypothesis']`
        field.
    *IS_DIRECT*
        Whether the Statement represents a direct physical interactions,
        as recorded by the `evidence.epistemics['direct']` field.

    In addition, if the `add_curation_cols` flag is set when calling
    :py:meth:`TsvAssembler.make_model`, the following additional (empty)
    columns will be added, to be filled out by curators:

    *AG_A_IDS_CORRECT*
        Correctness of Agent A grounding.
    *AG_A_STATE_CORRECT*
        Correctness of Agent A context (e.g., modification, bound, and other
        conditions).
    *AG_B_IDS_CORRECT, AG_B_STATE_CORRECT*
        As above, for Agent B.
    *EVENT_CORRECT*
        Whether the event is supported by the evidence text if the entities
        (Agents A and B) are considered as placeholders (i.e.,
        ignoring the correctness of their grounding).
    *RES_CORRECT*
        For Modification statements, whether the amino acid residue indicated
        by the Statement is supported by the evidence.
    *POS_CORRECT*
        For Modification statements, whether the amino acid position indicated
        by the Statement is supported by the evidence.
    *SUBJ_ACT_CORRECT*
        For Activation/Inhibition Statements, whether the activity indicated
        for the subject (Agent A) is supported by the evidence.
    *OBJ_ACT_CORRECT*
        For Activation/Inhibition Statements, whether the activity indicated
        for the object (Agent B) is supported by the evidence.
    *HYP_CORRECT*
        Whether the Statement is correctly flagged as a hypothesis.
    *HYP_CORRECT*
        Whether the Statement is correctly flagged as direct.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be assembled.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.
    """
    def __init__(self, statements=None):
        if not statements:
            self.statements = []
        else:
            self.statements = statements

    def add_statements(self, stmts):
        self.statements.extend(stmts)

    def make_model(self, output_file, add_curation_cols=False, up_only=False):
        """Export the statements into a tab-separated text file.

        Parameters
        ----------
        output_file : str
            Name of the output file.
        add_curation_cols : bool
            Whether to add columns to facilitate statement curation. Default
            is False (no additional columns).
        up_only : bool
            Whether to include identifiers.org links *only* for the Uniprot
            grounding of an agent when one is available. Because most
            spreadsheets allow only a single hyperlink per cell, this can makes
            it easier to link to Uniprot information pages for curation
            purposes. Default is False.
        """
        stmt_header = ['INDEX', 'UUID', 'TYPE', 'STR',
                       'AG_A_TEXT', 'AG_A_LINKS', 'AG_A_STR',
                       'AG_B_TEXT', 'AG_B_LINKS', 'AG_B_STR',
                       'PMID', 'TEXT', 'IS_HYP', 'IS_DIRECT']
        if add_curation_cols:
            stmt_header = stmt_header + \
                          ['AG_A_IDS_CORRECT', 'AG_A_STATE_CORRECT',
                           'AG_B_IDS_CORRECT', 'AG_B_STATE_CORRECT',
                           'EVENT_CORRECT',
                           'RES_CORRECT', 'POS_CORRECT', 'SUBJ_ACT_CORRECT',
                           'OBJ_ACT_CORRECT', 'HYP_CORRECT', 'DIRECT_CORRECT']
        rows = [stmt_header]

        for ix, stmt in enumerate(self.statements):
            # Complexes
            if len(stmt.agent_list()) > 2:
                logger.info("Skipping statement with more than two members: %s"
                            % stmt)
                continue
            # Self-modifications, ActiveForms
            elif len(stmt.agent_list()) == 1:
                ag_a = stmt.agent_list()[0]
                ag_b = None
            # All others
            else:
                (ag_a, ag_b) = stmt.agent_list()
            # Put together the data row
            row = [ix+1, stmt.uuid, stmt.__class__.__name__, str(stmt)] + \
                  _format_agent_entries(ag_a, up_only) + \
                  _format_agent_entries(ag_b, up_only) + \
                  [stmt.evidence[0].pmid, stmt.evidence[0].text,
                   stmt.evidence[0].epistemics.get('hypothesis', ''),
                   stmt.evidence[0].epistemics.get('direct', '')]
            if add_curation_cols:
                row = row + ([''] * 11)
            rows.append(row)
        # Write to file
        write_unicode_csv(output_file, rows, delimiter='\t')


def _format_id(ns, id):
    """Format a namespace/ID pair for display and curation."""
    label = '%s:%s' % (ns, id)
    label = label.replace(' ', '_')
    url = get_identifiers_url(ns, id)
    return (label, url)


def _format_agent_entries(agent, up_only):
    if agent is None:
        return ['', '', '']
    # Agent text/name
    agent_text = agent.db_refs.get('TEXT')
    if agent_text is None:
        agent_text = agent.name
    # Agent db_refs str
    db_refs = copy(agent.db_refs)
    if 'TEXT' in db_refs:
        db_refs.pop('TEXT')
    db_refs_str = ','.join(['%s|%s' % (k, v)
                            for k, v in db_refs.items()])
    # Agent links
    identifier_links = []
    if up_only and 'UP' in db_refs:
        up_label, up_url = _format_id('UP', db_refs['UP'])
        identifier_links = [up_url]
    else:
        for ns, id in db_refs.items():
            label, url = _format_id(ns, id)
            if url is None:
                identifier_links.append(label)
            else:
                identifier_links.append(url)
    links_str = ', '.join(identifier_links)
    return [agent_text, links_str, str(agent)]

