"""
Format a set of INDRA Statements into an HTML-formatted report which also
supports curation.
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import re
import uuid
import logging
import itertools
from collections import OrderedDict
from os.path import abspath, dirname, join, exists, getmtime, sep

from jinja2 import Environment, BaseLoader, TemplateNotFound, FileSystemLoader

from indra.statements import *
from indra.assemblers.english import EnglishAssembler
from indra.databases import get_identifiers_url
from indra.util.statement_presentation import group_and_sort_statements, \
    make_string_from_sort_key, make_top_level_label_from_names_key

logger = logging.getLogger(__name__)
HERE = dirname(abspath(__file__))


loader = FileSystemLoader(join(HERE, 'templates'))
env = Environment(loader=loader)

default_template = env.get_template('indra/statements_view.html')

color_schemes = {
    'dark': ['#b2df8a', '#000099', '#6a3d9a', '#1f78b4', '#fdbf6f', '#ff7f00',
             '#cab2d6', '#fb9a99', '#a6cee3', '#33a02c', '#b15928', '#e31a1c'],
    'light': ['#bc80bd', '#fccde5', '#b3de69', '#80b1d3', '#fb8072', '#bebada',
              '#fdb462', '#8dd3c7', '#ffffb3', '#d9d9d9', '#ccebc5', '#ffed6f']
}


def color_gen(scheme):
    while True:
        for color in color_schemes[scheme]:
            yield color


SOURCE_COLORS = [
    ('databases', {'color': 'black',
                   'sources': dict(zip(['phosphosite', 'cbn', 'pc11',
                                        'biopax', 'bel_lc',
                                        'signor', 'biogrid', 'tas',
                                        'lincs_drug', 'hprd', 'trrust'],
                                       color_gen('light')))}),
    ('reading', {'color': 'white',
                 'sources': dict(zip(['geneways', 'tees', 'isi', 'trips',
                                      'rlimsp', 'medscan', 'sparser', 'reach'],
                                     color_gen('light')))}),
]

SRC_KEY_DICT = {src: src for _, d in SOURCE_COLORS
                for src in d['sources'].keys()}


class HtmlAssembler(object):
    """Generates an HTML-formatted report from INDRA Statements.

    The HTML report format includes statements formatted in English
    (by the EnglishAssembler), text and metadata for the Evidence
    object associated with each Statement, and a Javascript-based curation
    interface linked to the INDRA database (access permitting). The interface
    allows for curation of statements at the evidence level by letting the
    user specify type of error and (optionally) provide a short description of
    of the error.

    Parameters
    ----------
    statements : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler. Statements
        can also be added using the add_statements method after the assembler
        has been instantiated.
    summary_metadata : Optional[dict]
        Dictionary of statement corpus metadata such as that provided by the
        INDRA REST API. Default is None. Each value should be a concise
        summary of O(1), not of order the length of the list, such as the
        evidence totals. The keys should be informative human-readable strings.
    ev_totals : Optional[dict]
        A dictionary of the total evidence available for each
        statement indexed by hash. Default: None
    source_counts : Optional[dict]
        A dictionary of the itemized evidence counts, by source, available for
        each statement, indexed by hash. Default: None.
    title : str
        The title to be printed at the top of the page.
    db_rest_url : Optional[str]
        The URL to a DB REST API to use for links out to further evidence.
        If given, this URL will be prepended to links that load additional
        evidence for a given Statement. One way to obtain this value is from
        the configuration entry indra.config.get_config('INDRA_DB_REST_URL').
        If None, the URLs are constructed as relative links.
        Default: None

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to assemble.
    model : str
        The HTML report formatted as a single string.
    metadata : dict
        Dictionary of statement list metadata such as that provided by the
        INDRA REST API.
    ev_totals : dict
        A dictionary of the total evidence available for each
        statement indexed by hash.
    db_rest_url : str
        The URL to a DB REST API.
    """

    def __init__(self, statements=None, summary_metadata=None, ev_totals=None,
                 source_counts=None, title='INDRA Results', db_rest_url=None):
        self.title = title
        self.statements = [] if statements is None else statements
        self.metadata = {} if summary_metadata is None \
            else summary_metadata
        self.ev_totals = {} if ev_totals is None else ev_totals
        self.source_counts = {} if source_counts is None else source_counts
        self.db_rest_url = db_rest_url
        self.model = None

    def add_statements(self, statements):
        """Add a list of Statements to the assembler.

        Parameters
        ----------
        statements : list[indra.statements.Statement]
            A list of INDRA Statements to be added to the assembler.
        """
        self.statements += statements

    def make_model(self, template=None, **template_kwargs):
        """Return the assembled HTML content as a string.

        Returns
        -------
        str
            The assembled HTML as a string.
        """
        # Get an iterator over the statements, carefully grouped.
        stmt_rows = group_and_sort_statements(
            self.statements,
            self.ev_totals if self.ev_totals else None,
            self.source_counts if self.source_counts else None)

        # Do some extra formatting.
        tl_stmts = OrderedDict()
        for row in stmt_rows:
            # Distinguish between the cases with
            if self.source_counts:
                key, verb, stmts, tl_counts, src_counts = row
            else:
                key, verb, stmts = row
                src_counts = None
                tl_counts = None

            names = key[1]
            tl_key = '-'.join([str(name) for name in names])
            tl_label = make_top_level_label_from_names_key(names)

            if tl_key not in tl_stmts.keys():
                tl_stmts[tl_key] = {'html_key': str(uuid.uuid4()),
                                    'label': tl_label,
                                    'source_counts': tl_counts,
                                    'stmts_formatted': []}

            # This will now be ordered by prevalence and entity pairs.
            stmt_info_list = []
            for stmt in stmts:
                stmt_hash = stmt.get_hash(shallow=True)
                ev_list = self._format_evidence_text(stmt)
                english = self._format_stmt_text(stmt)
                if self.ev_totals:
                    tot_ev = self.ev_totals.get(int(stmt_hash), '?')
                    if tot_ev == '?':
                        logger.warning('The hash %s was not found in the '
                                       'evidence totals dict.' % stmt_hash)
                    evidence_count_str = '%s / %s' % (len(ev_list), tot_ev)
                else:
                    evidence_count_str = str(len(ev_list))
                stmt_info_list.append({
                    'hash': stmt_hash,
                    'english': english,
                    'evidence': ev_list,
                    'evidence_count': evidence_count_str,
                    'source_count': self.source_counts.get(stmt_hash)})

            # Generate the short name for the statement and a unique key.
            short_name = make_string_from_sort_key(key, verb)
            short_name_key = str(uuid.uuid4())

            new_tpl = (short_name, short_name_key, stmt_info_list, src_counts)
            tl_stmts[tl_key]['stmts_formatted'].append(new_tpl)

        metadata = {k.replace('_', ' ').title(): v
                    for k, v in self.metadata.items()
                    if not isinstance(v, list) and not isinstance(v, dict)}
        if self.db_rest_url and not self.db_rest_url.endswith('statements'):
            db_rest_url = self.db_rest_url + '/statements'
        else:
            db_rest_url = None

        # Fill the template.
        if template is None:
            template = default_template
        self.model = template.render(stmt_data=tl_stmts,
                                     metadata=metadata, title=self.title,
                                     db_rest_url=db_rest_url,
                                     source_colors=SOURCE_COLORS,
                                     source_key_dict=SRC_KEY_DICT,
                                     **template_kwargs)
        return self.model

    def append_warning(self, msg):
        """Append a warning message to the model to expose issues."""
        assert self.model is not None, "You must already have run make_model!"
        addendum = ('\t<span style="color:red;">(CAUTION: %s occurred when '
                    'creating this page.)</span>' % msg)
        self.model = self.model.replace(self.title, self.title + addendum)
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
               isinstance(stmt, ActiveForm) or isinstance(stmt, Conversion) or\
               isinstance(stmt, Translocation):
                return 'other'
            else:
                assert len(stmt.agent_list()) == 2, (len(stmt.agent_list()),
                                                     type(stmt))
                return 'subject' if ag_ix == 0 else 'object'

        ev_list = []
        for ix, ev in enumerate(stmt.evidence):
            # Expand the source api to include the sub-database
            if ev.source_api == 'biopax' and \
               'source_sub_id' in ev.annotations and \
               ev.annotations['source_sub_id']:
               source_api = '%s:%s' % (ev.source_api,
                                       ev.annotations['source_sub_id'])
            else:
                source_api = ev.source_api
            # Prepare the evidence text
            if ev.text is None:
                format_text = None
            else:
                indices = []
                for ix, ag in enumerate(stmt.agent_list()):
                    if ag is None:
                        continue
                    # If the statement has been preassembled, it will have
                    # this entry in annotations
                    try:
                        ag_text = ev.annotations['agents']['raw_text'][ix]
                        if ag_text is None:
                            raise KeyError
                    # Otherwise we try to get the agent text from db_refs
                    except KeyError:
                        ag_text = ag.db_refs.get('TEXT')
                    if ag_text is None:
                        continue
                    role = get_role(ix)
                    # Get the tag with the correct badge
                    tag_start = '<span class="badge badge-%s">' % role
                    tag_close = '</span>'
                    # Build up a set of indices
                    indices += [(m.start(), m.start() + len(ag_text),
                                 ag_text, tag_start, tag_close)
                                 for m in re.finditer(re.escape(ag_text),
                                                      ev.text)]
                format_text = tag_text(ev.text, indices)

            ev_list.append({'source_api': source_api,
                            'pmid': ev.pmid,
                            'text_refs': ev.text_refs,
                            'text': format_text,
                            'source_hash': ev.source_hash })

        return ev_list

    @staticmethod
    def _format_stmt_text(stmt):
        # Get the English assembled statement
        ea = EnglishAssembler([stmt])
        english = ea.make_model()
        if not english:
            english = str(stmt)
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
            # FIXME: the EnglishAssembler capitalizes the first letter of
            # each sentence. In some cases this causes no match here
            # and not produce agent links.
            indices += [(m.start(), m.start() + len(ag.name), ag.name,
                         tag_start, tag_close)
                         for m in re.finditer(re.escape(ag.name), english)]
        return tag_text(english, indices)


def id_url(ag):
    # Return identifier URLs in a prioritized order
    for db_name in ('FPLX', 'HGNC', 'UP',
                    'GO', 'MESH',
                    'CHEBI', 'PUBCHEM', 'HMDB',
                    'IP', 'PF', 'NXPFA',
                    'MIRBASEM', 'MIRBASE',
                    'NCIT',
                    'UN', 'HUME', 'CWMS', 'SOFIA'):
        if db_name in ag.db_refs:
            # Handle a special case where a list of IDs is given
            if isinstance(ag.db_refs[db_name], list):
                db_id = ag.db_refs[db_name][0]
                if db_name == 'CHEBI':
                    if not db_id.startswith('CHEBI'):
                        db_id = 'CHEBI:%s' % db_id
                elif db_name in ('UN', 'HUME'):
                    db_id = db_id[0]
            else:
                db_id = ag.db_refs[db_name]
            return get_identifiers_url(db_name, db_id)


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
