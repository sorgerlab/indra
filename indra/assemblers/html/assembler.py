"""
Format a set of INDRA Statements into an HTML-formatted report which also
supports curation.
"""

import re
import uuid
import logging
import itertools
from collections import OrderedDict
from os.path import abspath, dirname, join

from jinja2 import Environment, FileSystemLoader

from indra.statements import *
from indra.sources import SOURCE_INFO
from indra.statements.agent import default_ns_order
from indra.statements.validate import validate_id
from indra.databases.identifiers import get_identifiers_url, ensure_prefix
from indra.assemblers.english import EnglishAssembler, AgentWithCoordinates
from indra.util.statement_presentation import group_and_sort_statements, \
    make_top_level_label_from_names_key, make_stmt_from_sort_key, \
    reader_sources, db_sources, all_sources, get_available_source_counts, \
    get_available_ev_counts, standardize_counts
from indra.literature import id_lookup

logger = logging.getLogger(__name__)
HERE = dirname(abspath(__file__))


loader = FileSystemLoader(join(HERE, 'templates'))
env = Environment(loader=loader)

default_template = env.get_template('indra/statements_view.html')

color_schemes = {
    'dark': ['#b2df8a', '#000099', '#6a3d9a', '#1f78b4', '#fdbf6f', '#ff7f00',
             '#cab2d6', '#fb9a99', '#a6cee3', '#33a02c', '#b15928', '#e31a1c'],
    'light': ['#bebada', '#fdb462', '#b3de69', '#80b1d3', '#bc80bd', '#fccde5',
              '#fb8072', '#d9d9d9', '#8dd3c7', '#ffed6f', '#ccebc5', '#e0e03d',
              '#ffe8f4', '#acfcfc', '#dd99ff', '#00d4a6']
}


def color_gen(scheme):
    while True:
        for color in color_schemes[scheme]:
            yield color


def make_source_colors(databases, readers):
    rdr_ord = ['reach', 'sparser', 'medscan', 'trips', 'eidos']
    readers.sort(key=lambda r: rdr_ord.index(r) if r in rdr_ord else len(rdr_ord))
    reader_colors_list = list(zip(readers, color_gen('light')))
    reader_colors_list.reverse()
    reader_colors = dict(reader_colors_list)
    db_colors = dict(zip(databases, color_gen('light')))
    return [('databases', {'color': 'black', 'sources': db_colors}),
            ('reading', {'color': 'white', 'sources': reader_colors})]


DEFAULT_SOURCE_COLORS = make_source_colors(db_sources, reader_sources)


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
    ev_counts : Optional[dict]
        A dictionary of the total evidence available for each
        statement indexed by hash. If not provided, the statements that are
        passed to the constructor are used to determine these, with whatever
        evidences these statements carry.
    ev_totals : Optional[dict]
        DEPRECATED. Same as ev_counts which should be used instead.
    source_counts : Optional[dict]
        A dictionary of the itemized evidence counts, by source, available for
        each statement, indexed by hash. If not provided, the statements
        that are passed to the constructor are used to determine these, with
        whatever evidences these statements carry.
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
    ev_counts : dict
        A dictionary of the total evidence available for each
        statement indexed by hash.
    db_rest_url : str
        The URL to a DB REST API.
    """

    def __init__(self, statements=None, summary_metadata=None, ev_totals=None,
                 ev_counts=None, source_counts=None, curation_dict=None,
                 title='INDRA Results', db_rest_url=None):
        self.title = title
        self.statements = [] if statements is None else statements
        self.metadata = {} if summary_metadata is None \
            else summary_metadata
        # If the deprecated parameter is used, we make sure we take it
        if not ev_counts and ev_totals:
            ev_counts = ev_totals
        self.ev_counts = get_available_ev_counts(self.statements) \
            if ev_counts is None else standardize_counts(ev_counts)
        self.source_counts = get_available_source_counts(self.statements) \
            if source_counts is None else standardize_counts(source_counts)
        self.curation_dict = {} if curation_dict is None else curation_dict
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

    def make_json_model(self, with_grouping=True, no_redundancy=False):
        """Return the JSON used to create the HTML display.

        Parameters
        ----------
        with_grouping : Optional[bool]
            If True, statements will be grouped under multiple sub-headings. If
            False, all headings will be collapsed into one on every level, with
            all statements placed under a single heading. Default: False
        no_redundancy : Optional[bool]
            If True, any group of statements that was already presented under
            a previous heading will be skipped. This is typically the case
            for complexes where different permutations of complex members
            are presented. By setting this argument to True, these can be
            eliminated. Default: False

        Returns
        -------
        json : dict
            A complexly structured JSON dict containing grouped statements and
            various metadata.
        """
        # Get an iterator over the statements, carefully grouped.
        stmt_rows = group_and_sort_statements(
            self.statements,
            self.ev_counts if self.ev_counts else None,
            self.source_counts if self.source_counts else None)

        # Do some extra formatting.
        stmts = OrderedDict()
        agents = {}
        previous_stmt_set = set()
        all_previous_stmts = set()
        for row in stmt_rows:
            # Distinguish between the cases with source counts and without.
            if self.source_counts:
                key, verb, stmts_group, tl_counts, src_counts = row
            else:
                key, verb, stmts_group = row
                src_counts = None
                tl_counts = None
            curr_stmt_set = {s.get_hash() for s in stmts_group}
            if curr_stmt_set == previous_stmt_set:
                continue
            elif no_redundancy and curr_stmt_set <= all_previous_stmts:
                continue
            else:
                previous_stmt_set = curr_stmt_set
                all_previous_stmts |= curr_stmt_set

            # We will keep track of some of the meta data for this stmt group.
            # NOTE: Much of the code relies heavily on the fact that the Agent
            # objects in `meta_agents` are references to the Agent's in the
            # Statement object `meta_stmts`.
            meta_agents = []
            meta_stmt = make_stmt_from_sort_key(key, verb, meta_agents)
            meta_agent_dict = {ag.name: ag for ag in meta_agents
                               if ag is not None}

            # This will now be ordered by prevalence and entity pairs.
            stmt_info_list = []
            for stmt in stmts_group:
                stmt_hash = stmt.get_hash(shallow=True)

                # Try to accumulate db refs in the meta agents.
                for ag in stmt.agent_list():
                    if ag is None:
                        continue
                    # Get the corresponding meta-agent
                    meta_ag = meta_agent_dict.get(ag.name)
                    if not meta_ag:
                        continue
                    _cautiously_merge_refs(ag, meta_ag)

                # Format some strings nicely.
                ev_list = _format_evidence_text(stmt, self.curation_dict)
                english = _format_stmt_text(stmt)
                if self.ev_counts:
                    tot_ev = self.ev_counts.get(int(stmt_hash), '?')
                    if tot_ev == '?':
                        logger.warning('The hash %s was not found in the '
                                       'evidence totals dict.' % stmt_hash)
                    evidence_count_str = '%s / %s' % (len(ev_list), tot_ev)
                else:
                    evidence_count_str = str(len(ev_list))

                stmt_info_list.append({
                    'hash': str(stmt_hash),
                    'english': english,
                    'evidence': ev_list,
                    'evidence_count': evidence_count_str,
                    'source_count': self.source_counts.get(stmt_hash)})

            # Clean out invalid fields from the meta agents.
            for ag in meta_agents:
                if ag is None:
                    continue
                for dbn, dbid in list(ag.db_refs.items()):
                    if isinstance(dbid, set):
                        logger.info("Removing %s from refs due to too many "
                                    "matches: %s" % (dbn, dbid))
                        del ag.db_refs[dbn]

            # Update the top level grouping.
            if isinstance(stmt, (ActiveForm, HasActivity)):
                tl_names = [key[1][0]]
            elif isinstance(stmt, Conversion):
                tl_names = [key[1][0]] + [*key[1][1]] + [*key[1][2]]
            else:
                tl_names = key[1]
            if with_grouping:
                tl_key = '-'.join([str(name) for name in tl_names])
                tl_agents = {name: Agent(name) for name in tl_names
                             if name is not None}
                for ag in tl_agents.values():
                    meta_ag = meta_agent_dict.get(ag.name)
                    if meta_ag is None:
                        continue
                    ag.db_refs.update(meta_ag.db_refs)
                tl_label = None
            else:
                tl_key = 'all-statements'
                tl_label = 'All Statements'
                tl_agents = None

            if tl_key not in stmts.keys():
                agents[tl_key] = tl_agents
                stmts[tl_key] = {'html_key': str(uuid.uuid4()),
                                 'source_counts': tl_counts,
                                 'stmts_formatted': [],
                                 'names': tl_names}
                if tl_label:
                    stmts[tl_key]['label'] = tl_label
            elif with_grouping:
                for name, existing_ag in agents[tl_key].items():
                    new_ag = tl_agents.get(name)
                    if new_ag is None:
                        continue
                    _cautiously_merge_refs(new_ag, existing_ag)

            # Generate the short name for the statement and a unique key.
            existing_list = stmts[tl_key]['stmts_formatted']
            if with_grouping or not existing_list:
                if with_grouping:
                    # See note above: this is where the work on meta_agents is
                    # applied because the agents are references.
                    short_name = _format_stmt_text(meta_stmt)
                    short_name_key = str(uuid.uuid4())
                else:
                    short_name = "All Statements Sub Group"
                    short_name_key = "all-statements-sub-group"
                new_dict = {'short_name': short_name,
                            'short_name_key': short_name_key,
                            'stmt_info_list': stmt_info_list,
                            'src_counts': src_counts}
                existing_list.append(new_dict)
            else:
                existing_list[0]['stmt_info_list'].extend(stmt_info_list)
                if src_counts:
                    existing_list[0]['src_counts'].update(src_counts)

        # Add labels for each top level group (tlg).
        if with_grouping:
            for tl_key, tlg in stmts.items():
                tl_agents = list(agents[tl_key].values())
                for ag in tl_agents:
                    for dbn, dbid in list(ag.db_refs.items()):
                        if isinstance(dbid, set):
                            logger.info("Removing %s from top level refs "
                                        "due to multiple matches: %s"
                                        % (dbn, dbid))
                            del ag.db_refs[dbn]
                tl_label = make_top_level_label_from_names_key(tlg['names'])
                tl_label = re.sub("<b>(.*?)</b>", r"\1", tl_label)
                tl_label = tag_agents(tl_label, tl_agents)
                tlg['label'] = tl_label

        return stmts

    def make_model(self, template=None, with_grouping=True,
                   add_full_text_search_link=False, no_redundancy=False,
                   **template_kwargs):
        """Return the assembled HTML content as a string.

        Parameters
        ----------
        template : a Template object
            Manually pass a Jinja template to be used in generating the HTML.
            The template is responsible for rendering essentially the output of
            `make_json_model`.
        with_grouping : bool
            If True, statements will be grouped under multiple sub-headings. If
            False, all headings will be collapsed into one on every level, with
            all statements placed under a single heading.
        add_full_text_search_link : bool
            If True, link with Text fragment search in PMC journal will be
            added for the statements.  
        no_redundancy : Optional[bool]
            If True, any group of statements that was already presented under
            a previous heading will be skipped. This is typically the case
            for complexes where different permutations of complex members
            are presented. By setting this argument to True, these can be
            eliminated. Default: False

            All other keyword arguments are passed along to the template. If you
            are using a custom template with args that are not passed below, this
            is how you pass them.

        Returns
        -------
        str
            The assembled HTML as a string.
        """
        tl_stmts = self.make_json_model(with_grouping,
                                        no_redundancy=no_redundancy)

        if add_full_text_search_link:
            for statement in tl_stmts:
                statement = tl_stmts[statement]
                for stmt_formatted in statement["stmts_formatted"]:
                    for stmt_info in stmt_formatted["stmt_info_list"]:
                        for evidence in stmt_info["evidence"]:
                            if 'PMCID' not in evidence.get('text_refs', {}):
                                if evidence.get('pmid'):
                                    ev_pmcid = id_lookup(
                                        evidence['pmid'], 'pmid') \
                                        .get('pmcid', None)
                                    if ev_pmcid:
                                        evidence['pmcid'] = ev_pmcid
                            else:
                                evidence['pmcid'] = \
                                    evidence['text_refs']['PMCID']

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
        if self.source_counts and 'source_key_dict' not in template_kwargs:
            template_kwargs['source_key_dict'] = \
                {src: src for src in all_sources}
        if 'source_colors' not in template_kwargs:
            template_kwargs['source_colors'] = DEFAULT_SOURCE_COLORS
        if 'source_info' not in template_kwargs:
            template_kwargs['source_info'] = SOURCE_INFO.copy()
        if 'simple' not in template_kwargs:
            template_kwargs['simple'] = True

        self.model = template.render(stmt_data=tl_stmts,
                                     metadata=metadata, title=self.title,
                                     db_rest_url=db_rest_url,
                                     add_full_text_search_link=add_full_text_search_link,  # noqa
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


def _format_evidence_text(stmt, curation_dict=None, correct_tags=None):
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
    if curation_dict is None:
        curation_dict = {}
    if correct_tags is None:
        correct_tags = ['correct']

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
                            for m in re.finditer(re.escape(ag_text), ev.text)]
            format_text = tag_text(ev.text, indices)

        curation_key = (stmt.get_hash(), ev.source_hash)
        curations = curation_dict.get(curation_key, [])
        num_curations = len(curations)
        num_correct = len(
            [cur for cur in curations if cur['error_type'] in correct_tags])
        num_incorrect = num_curations - num_correct
        text_refs = {k.upper(): v for k, v in ev.text_refs.items()}
        ev_list.append({'source_api': source_api,
                        'pmid': ev.pmid,
                        'text_refs': text_refs,
                        'text': format_text,
                        'source_hash': str(ev.source_hash),
                        'num_curations': num_curations,
                        'num_correct': num_correct,
                        'num_incorrect': num_incorrect
                        })

    return ev_list


def _format_stmt_text(stmt):
    # Get the English assembled statement
    ea = EnglishAssembler([stmt])
    english = ea.make_model()
    if not english:
        english = str(stmt)
        return tag_agents(english, stmt.agent_list())
    return tag_agents(english, ea.stmt_agents[0])


def _cautiously_merge_refs(from_ag, to_ag):
    # Check the db refs for this agent against the meta agent
    for dbn, dbid in from_ag.db_refs.items():
        if dbn == 'TEXT':
            continue
        meta_dbid = to_ag.db_refs.get(dbn)
        if isinstance(meta_dbid, set):
            # If we've already marked this one add to the set.
            to_ag.db_refs[dbn].add(dbid)
        elif meta_dbid is not None and meta_dbid != dbid:
            # If we've seen it before and don't agree, mark it.
            to_ag.db_refs[dbn] = {to_ag.db_refs[dbn], dbid}
        elif meta_dbid is None:
            # Otherwise, add it.
            to_ag.db_refs[dbn] = dbid


def tag_agents(english, agents):
    # Agents can be AgentWithCoordinates (preferred) or regular Agent objects
    indices = []
    for ag in agents:
        if ag is None or not ag.name:
            continue
        url = id_url(ag)
        if url is None:
            tag_start = '<b>'
            tag_close = '</b>'
        else:
            tag_start = "<a href='%s' target='_blank'>" % url
            tag_close = "</a>"
        # If coordinates are passed, use them. Otherwise, try to find agent
        # names in english text
        if isinstance(ag, AgentWithCoordinates):
            index = (ag.coords[0], ag.coords[1], ag.name, tag_start, tag_close)
            indices.append(index)
        elif isinstance(ag, Agent):
            found = False
            for m in re.finditer(re.escape(ag.name), english):
                index = (m.start(), m.start() + len(ag.name), ag.name,
                         tag_start, tag_close)
                indices.append(index)
                found = True
            if not found and \
                    english.startswith(re.escape(ag.name).capitalize()):
                index = (0, len(ag.name), ag.name, tag_start, tag_close)
                indices.append(index)
    return tag_text(english, indices)


link_namespace_order = default_ns_order + \
    ['CHEMBL', 'DRUGBANK', 'PUBCHEM', 'HMDB', 'HMS-LINCS', 'CAS',
     'IP', 'PF', 'NXPFA', 'MIRBASEM', 'NCIT', 'WM']


def id_url(ag):
    # Return identifier URLs in a prioritized order
    # TODO: we should add handling for UPPRO here, however, that would require
    # access to UniProt client resources in the context of the DB REST API
    # which could be problematic
    for db_name in link_namespace_order:
        if db_name in ag.db_refs:
            # Handle a special case where a list of IDs is given
            if isinstance(ag.db_refs[db_name], list):
                db_id = ag.db_refs[db_name][0]
                if db_name == 'WM':
                    db_id = db_id[0]
            else:
                db_id = ag.db_refs[db_name]
            # We can add more name spaces here if there are issues
            if db_name in {'CHEBI'}:
                db_id = ensure_prefix('CHEBI', db_id)
            # Here we validate IDs to make sure we don't surface invalid
            # links.
            if not validate_id(db_name, db_id):
                logger.debug('Invalid grounding encountered: %s:%s' %
                             (db_name, db_id))
                continue
            # Finally, we return a valid identifiers.org URL
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
        # Capitalize if it's the beginning of a sentence
        if i == 0:
            ag_text = ag_text[0].upper() + ag_text[1:]
        # Add the text before this agent, if any
        format_text += text[start_pos:i]
        # Add wrapper for this entity
        format_text += tag_start + ag_text + tag_close
        # Now set the next start position
        start_pos = j
    # Add the last section of text
    format_text += text[start_pos:]
    return format_text
