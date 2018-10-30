"""
Format a set of INDRA Statements into an HTML-based format.
"""
import sys
import csv
from os.path import abspath, dirname, join
import re
import json
from collections import defaultdict
from jinja2 import Template
from indra.sources import indra_db_rest
from indra.assemblers.english import EnglishAssembler
from indra.statements import *
from indra.databases import get_identifiers_url

# Create a template object from the template file, load once
template_path = join(dirname(abspath(__file__)), 'template.html')
with open(template_path, 'rt') as f:
    template_str = f.read()
    template = Template(template_str)

# TODO:
# - Highlight text in english assembled sentences
# - For both, add links to identifiers.org


class HtmlAssembler(object):
    """Generates an HTML-formatted report from INDRA Statements.

    The HTML report format includes statements formatted in English
    (formatted by the EnglishAssembler), text and metadata for the Evidence
    object associated with each Statement, and a Javascript-based curation
    interface linked to the INDRA database (access permitting).

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler.
    rest_api_results : dict
        Dictionary of query metadata provided by the INDRA REST API. Default
        is None.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to assemble.
    model : str
        The HTML report formatted as a single string.
    """
    def __init__(self, stmts=None, rest_api_results=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        if rest_api_results is None:
            self.rest_api_results = {}
        else:
            self.rest_api_results = rest_api_results
        self.model = None

    def make_model(self):
        stmts_formatted = []
        for stmt in self.statements:
            stmt_hash = stmt.get_hash(shallow=True)
            ea = EnglishAssembler([stmt])
            english = ea.make_model()
            ev_list = self.format_evidence_text(stmt)
            if self.rest_api_results:
                total_evidence = self.rest_api_results['evidence_totals']\
                                                      [int(stmt_hash)]
                evidence_count_str = '%s / %s' % (len(ev_list), total_evidence)
            else:
                evidence_count_str = str(len(ev_list))
            stmts_formatted.append({
                'hash': stmt_hash,
                'english': english,
                'evidence': ev_list,
                'evidence_count': evidence_count_str})
        return template.render(statements=stmts_formatted,
                               rest_api_results=self.rest_api_results)


    def format_evidence_text(self, stmt):
        """Highlight subject and object in raw text strings."""
        def get_role(ag_ix):
            if isinstance(stmt, Complex) or \
               isinstance(stmt, SelfModification) or \
               isinstance(stmt, ActiveForm) or isinstance(stmt, Conversion) or \
               isinstance(stmt, Translocation):
                return 'other'
            else:
                assert len(stmt.agent_list()) == 2
                return 'subject' if ag_ix == 0 else 'object'

        ev_list = []
        for ix, ev in enumerate(stmt.evidence):
            if ev.text is None:
                format_text = '(None available)'
            else:
                indices = []
                for ix, ag in enumerate(stmt.agent_list()):
                    ag_text = ev.annotations['agents']['raw_text'][ix]
                    if ag_text is None:
                        continue
                    # Build up a set of indices
                    indices += [(m.start(), m.start() + len(ag_text),
                                 ag_text, ix)
                                 for m in re.finditer(re.escape(ag_text),
                                                      ev.text)]
                # Sort the indices by their start position
                indices.sort(key=lambda x: x[0])
                # Now, add the marker text for each occurrence of the strings
                format_text = ''
                start_pos = 0
                for i, j, ag_text, ag_ix in indices:
                    role = get_role(ag_ix)
                    # Get the tag with the correct badge
                    tag_start = '<span class="label label-%s">' % role
                    tag_close = '</span>'
                    # Add the text before this agent, if any
                    format_text += ev.text[start_pos:i]
                    # Add wrapper for this entity
                    format_text += tag_start + ag_text + tag_close
                    # Now set the next start position
                    start_pos = j
                # Add the last section of text
                format_text += ev.text[start_pos:]
            ev_list.append({'source_api': ev.source_api,
                            'pmid': ev.pmid, 'text': format_text })
        return ev_list


