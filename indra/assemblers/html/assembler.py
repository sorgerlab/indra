"""
Format a set of INDRA Statements into an HTML-formatted report which also
supports curation.
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import itertools
from os.path import abspath, dirname, join
from jinja2 import Template
from indra.statements import *
from indra.assemblers.english import EnglishAssembler
from indra.databases import get_identifiers_url


# Create a template object from the template file, load once
template_path = join(dirname(abspath(__file__)), 'template.html')
with open(template_path, 'rt') as f:
    template_str = f.read()
    template = Template(template_str)


class HtmlAssembler(object):
    """Generates an HTML-formatted report from INDRA Statements.

    The HTML report format includes statements formatted in English
    (by the EnglishAssembler), text and metadata for the Evidence
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
    rest_api_results : dict
        Dictionary of query metadata provided by the INDRA REST API.
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
        """Return the assembled HTML content as a string.

        Returns
        -------
        str
            The assembled HTML as a string.
        """
        stmts_formatted = []
        for stmt in self.statements:
            stmt_hash = stmt.get_hash(shallow=True)
            ev_list = self._format_evidence_text(stmt)
            english = self._format_stmt_text(stmt)
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
        self.model = template.render(statements=stmts_formatted,
                                     rest_api_results=self.rest_api_results)
        return self.model

    def save_model(self, fname):
        """Save the assembled HTML into a file.

        Parameters
        ----------
        fname : str
            The path to the file to save the HTML into.
        """
        if self.model is None:
            self.make_model()

        with open(fname, 'wb') as fh:
            fh.write(self.model.encode('utf-8'))

    @staticmethod
    def _format_evidence_text(stmt):
        """Returns evidence metadata with highlighted evidence text.

        Parameters
        ----------
        stmt : indra.Statement
            The Statement with Evidence to be formatted.

        Returns
        -------
        list of dicts
            List of dictionaries corresponding to each Evidence object in the
            Statement's evidence list. Each dictionary has keys 'source_api',
            'pmid' and 'text', drawn from the corresponding fields in the
            Evidence objects. The text entry of the dict includes
            `<span>` tags identifying the agents referenced by the Statement.
        """
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
                    # If the statement has been preassembled, it will have
                    # this entry in annotations
                    try:
                        ag_text = ev.annotations['agents']['raw_text'][ix]
                    # Otherwise we try to get the agent text from db_refs
                    except KeyError:
                        ag_text = ag.db_refs.get('TEXT')
                    if ag_text is None:
                        continue
                    role = get_role(ix)
                    # Get the tag with the correct badge
                    tag_start = '<span class="label label-%s">' % role
                    tag_close = '</span>'
                    # Build up a set of indices
                    indices += [(m.start(), m.start() + len(ag_text),
                                 ag_text, tag_start, tag_close)
                                 for m in re.finditer(re.escape(ag_text),
                                                      ev.text)]
                format_text = tag_text(ev.text, indices)
            ev_list.append({'source_api': ev.source_api,
                            'text_refs': ev.text_refs,
                            'text': format_text,
                            'source_hash': ev.source_hash })

        return ev_list

    @staticmethod
    def _format_stmt_text(stmt):
        # Get the English assembled statement
        ea = EnglishAssembler([stmt])
        english = ea.make_model()
        indices = []
        for ag in stmt.agent_list():
            if ag is None or not ag.name:
                continue
            url = id_url(ag)
            if url is None:
                continue
            # Build up a set of indices
            tag_start = "<a href='%s'>" % url
            tag_close = "</a>"
            indices += [(m.start(), m.start() + len(ag.name), ag.name,
                         tag_start, tag_close)
                         for m in re.finditer(re.escape(ag.name), english)]
        return tag_text(english, indices)


def id_url(ag):
    # Return identifier URLs in a prioritized order
    for db_name in ('HGNC', 'FPLX', 'UP', 'IP', 'PF', 'NXPFA',
                    'MIRBASEM', 'MIRBASE',
                    'MESH', 'GO',
                    'HMDB', 'PUBCHEM', 'CHEBI',
                    'NCIT'):
        if db_name in ag.db_refs:
            return get_identifiers_url(db_name, ag.db_refs[db_name])


def tag_text(text, tag_info_list):
    """Apply start/end tags to spans of the given text.


    Parameters
    ----------
    text : str
        Text to be tagged
    tag_info_list : list of tuples
        Each tuple refers to a span of the given text. Fields are `(start_ix,
        end_ix, substring, start_tag, close_tag)`, where substring, start_tag,
        and close_tag are strings. If any of the given spans of text overlap,
        the longest span is used.

    Returns
    -------
    str
        String where the specified substrings have been surrounded by the
        given start and close tags.
    """
    # Check to tags for overlap and if there is any, return the subsumed
    # range. Return None if no overlap.
    def overlap(t1, t2):
        if range(max(t1[0], t2[0]), min(t1[1]-1, t2[1]-1)+1):
            if t1[1] - t1[0] >= t2[1] - t2[0]:
                return t2
            else:
                return t1
        else:
            return None
    # Remove subsumed tags
    for t1, t2 in list(itertools.combinations(tag_info_list, 2)):
        subsumed_tag = overlap(t1, t2)
        if subsumed_tag is not None:
            # Delete the subsumed tag from the list
            try:
                tag_ix = tag_info_list.index(subsumed_tag)
                del tag_info_list[tag_ix]
            # Ignore case where tag has already been deleted
            except ValueError:
                pass
    # Sort the indices by their start position
    tag_info_list.sort(key=lambda x: x[0])
    # Now, add the marker text for each occurrence of the strings
    format_text = ''
    start_pos = 0
    for i, j, ag_text, tag_start, tag_close in tag_info_list:
        # Add the text before this agent, if any
        format_text += text[start_pos:i]
        # Add wrapper for this entity
        format_text += tag_start + ag_text + tag_close
        # Now set the next start position
        start_pos = j
    # Add the last section of text
    format_text += text[start_pos:]
    return format_text


