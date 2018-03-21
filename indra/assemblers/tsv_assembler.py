from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from copy import copy
from indra.databases import get_identifiers_url
from indra.statements import *
from indra.util import write_unicode_csv
from indra.tools import assemble_corpus as ac

logger = logging.getLogger('tsv_assembler')

class TsvAssembler(object):
    """Assembles Statements into a set of tabular files for export or curation.

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

    def make_model(self, output_file, add_curation_cols=False):
        """Export the statements into a tab-separated text file.

        Parameters
        ----------
        output_file : str
            Name of the output file.
        add_curation_cols : bool
            Whether to add columns to facilitate statement curation. Default
            is False (no additional columns).
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
                  _format_agent_entries(ag_a) + _format_agent_entries(ag_b) + \
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


def _format_agent_entries(agent):
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
    if 'UP' in db_refs:
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

