__all__ = ['stmts_from_json', 'stmts_from_json_file', 'stmts_to_json',
           'stmts_to_json_file', 'draw_stmt_graph', 'pretty_print_stmts',
           'UnresolvedUuidError', 'InputError',
           'set_pretty_print_max_width']

import json
import logging
from typing import List, Optional

from indra.statements.statements import Statement, Unresolved


logger = logging.getLogger(__name__)


def stmts_from_json(json_in, on_missing_support='handle'):
    """Get a list of Statements from Statement jsons.

    In the case of pre-assembled Statements which have `supports` and
    `supported_by` lists, the uuids will be replaced with references to
    Statement objects from the json, where possible. The method of handling
    missing support is controled by the `on_missing_support` key-word argument.

    Parameters
    ----------
    json_in : iterable[dict]
        A json list containing json dict representations of INDRA Statements,
        as produced by the `to_json` methods of subclasses of Statement, or
        equivalently by `stmts_to_json`.
    on_missing_support : Optional[str]
        Handles the behavior when a uuid reference in `supports` or
        `supported_by` attribute cannot be resolved. This happens because uuids
        can only be linked to Statements contained in the `json_in` list, and
        some may be missing if only some of all the Statements from pre-
        assembly are contained in the list.

        Options:

        - *'handle'* : (default) convert unresolved uuids into `Unresolved`
          Statement objects.
        - *'ignore'* : Simply omit any uuids that cannot be linked to any
          Statements in the list.
        - *'error'* : Raise an error upon hitting an un-linkable uuid.

    Returns
    -------
    stmts : list[:py:class:`Statement`]
        A list of INDRA Statements.
    """

    stmts = []
    uuid_dict = {}
    for json_stmt in json_in:
        try:
            st = Statement._from_json(json_stmt)
        except Exception as e:
            logger.warning("Error creating statement: %s" % e)
            continue
        stmts.append(st)
        uuid_dict[st.uuid] = st
    for st in stmts:
        _promote_support(st.supports, uuid_dict, on_missing_support)
        _promote_support(st.supported_by, uuid_dict, on_missing_support)
    return stmts


def stmts_from_json_file(fname, format='json'):
    """Return a list of statements loaded from a JSON file.

    Parameters
    ----------
    fname : str
        Path to the JSON file to load statements from.
    format : Optional[str]
        One of 'json' to assume regular JSON formatting or
        'jsonl' assuming each statement is on a new line.

    Returns
    -------
    list[indra.statements.Statement]
        The list of INDRA Statements loaded from the JSOn file.
    """
    with open(fname, 'r') as fh:
        if format == 'json':
            return stmts_from_json(json.load(fh))
        else:
            return stmts_from_json([json.loads(line)
                                    for line in fh.readlines()])


def stmts_to_json_file(stmts, fname, format='json', **kwargs):
    """Serialize a list of INDRA Statements into a JSON file.

    Parameters
    ----------
    stmts : list[indra.statement.Statements]
        The list of INDRA Statements to serialize into the JSON file.
    fname : str
        Path to the JSON file to serialize Statements into.
    format : Optional[str]
        One of 'json' to use regular JSON with indent=1 formatting or
        'jsonl' to put each statement on a new line without indents.
    """
    sj = stmts_to_json(stmts, **kwargs)
    with open(fname, 'w') as fh:
        if format == 'json':
            json.dump(sj, fh, indent=1)
        else:
            for json_stmt in sj:
                json.dump(json_stmt, fh)
                fh.write('\n')


def stmts_to_json(stmts_in, use_sbo=False, matches_fun=None):
    """Return the JSON-serialized form of one or more INDRA Statements.

    Parameters
    ----------
    stmts_in : Statement or list[Statement]
        A Statement or list of Statement objects to serialize into JSON.
    use_sbo : Optional[bool]
        If True, SBO annotations are added to each applicable element of the
        JSON. Default: False
    matches_fun : Optional[function]
        A custom function which, if provided, is used to construct the
        matches key which is then hashed and put into the return value.
        Default: None

    Returns
    -------
    json_dict : dict
        JSON-serialized INDRA Statements.
    """
    if not isinstance(stmts_in, list):
        json_dict = stmts_in.to_json(use_sbo=use_sbo)
        return json_dict
    else:
        json_dict = [st.to_json(use_sbo=use_sbo, matches_fun=matches_fun)
                     for st in stmts_in]
    return json_dict


def _promote_support(sup_list, uuid_dict, on_missing='handle'):
    """Promote the list of support-related uuids to Statements, if possible."""
    valid_handling_choices = ['handle', 'error', 'ignore']
    if on_missing not in valid_handling_choices:
        raise InputError('Invalid option for `on_missing_support`: \'%s\'\n'
                         'Choices are: %s.'
                         % (on_missing, str(valid_handling_choices)))
    for idx, uuid in enumerate(sup_list):
        if uuid in uuid_dict.keys():
            sup_list[idx] = uuid_dict[uuid]
        elif on_missing == 'handle':
            sup_list[idx] = Unresolved(uuid)
        elif on_missing == 'ignore':
            sup_list.remove(uuid)
        elif on_missing == 'error':
            raise UnresolvedUuidError("Uuid %s not found in stmt jsons."
                                      % uuid)
    return


def draw_stmt_graph(stmts):
    """Render the attributes of a list of Statements as directed graphs.

    The layout works well for a single Statement or a few Statements at a time.
    This function displays the plot of the graph using plt.show().

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of one or more INDRA Statements whose attribute graph should
        be drawn.
    """
    import networkx
    try:
        import matplotlib.pyplot as plt
    except Exception:
        logger.error('Could not import matplotlib, not drawing graph.')
        return
    try:  # This checks whether networkx has this package to work with.
        import pygraphviz
    except Exception:
        logger.error('Could not import pygraphviz, not drawing graph.')
        return
    import numpy
    g = networkx.compose_all([stmt.to_graph() for stmt in stmts])
    plt.figure()
    plt.ion()
    g.graph['graph'] = {'rankdir': 'LR'}
    pos = networkx.drawing.nx_agraph.graphviz_layout(g, prog='dot')
    g = g.to_undirected()

    # Draw nodes
    options = {
        'marker': 'o',
        's': 200,
        'c': [0.85, 0.85, 1],
        'facecolor': '0.5',
        'lw': 0,
    }
    ax = plt.gca()
    nodelist = list(g)
    xy = numpy.asarray([pos[v] for v in nodelist])
    node_collection = ax.scatter(xy[:, 0], xy[:, 1], **options)
    node_collection.set_zorder(2)
    # Draw edges
    networkx.draw_networkx_edges(g, pos, arrows=False, edge_color='0.5')
    # Draw labels
    edge_labels = {(e[0], e[1]): e[2].get('label') for e in g.edges(data=True)}
    networkx.draw_networkx_edge_labels(g, pos, edge_labels=edge_labels)
    node_labels = {n[0]: n[1].get('label') for n in g.nodes(data=True)}
    for key, label in node_labels.items():
        if len(label) > 25:
            parts = label.split(' ')
            parts.insert(int(len(parts)/2), '\n')
            label = ' '.join(parts)
            node_labels[key] = label
    networkx.draw_networkx_labels(g, pos, labels=node_labels)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.show()


pretty_print_max_width = None


def set_pretty_print_max_width(new_max):
    """Set the max display width for pretty prints, in characters."""
    global pretty_print_max_width
    if new_max is not None and not isinstance(new_max, int):
        raise ValueError("Max width must be an integer or None.")
    pretty_print_max_width = new_max


def pretty_print_stmts(stmt_list: List[Statement],
                       stmt_limit: Optional[int] = None,
                       ev_limit: Optional[int] = 5,
                       width: Optional[int] = None) -> None:
    """Print a formatted list of statements along with evidence text.

    Requires the tabulate package (https://pypi.org/project/tabulate).

    Parameters
    ----------
    stmt_list : List[Statement]
        The list of INDRA Statements to be printed.
    stmt_limit : Optional[int]
        The maximum number of INDRA Statements to be printed. If None, all
        Statements are printed. (Default is None)
    ev_limit : Optional[int]
        The maximum number of Evidence to print for each Statement. If None, all
        evidence will be printed for each Statement. (Default is 5)
    width : Optional[int]
        Manually set the width of the table. If `None` the function will try to
        match the current terminal width using `os.get_terminal_size()`.  If
        this fails the width defaults to 80 characters. The maximum width can
        be controlled by setting :data:`pretty_print_max_width` using the
        :func:`set_pretty_print_max_width` function. This is useful in 
        Jupyter notebooks where the environment returns a terminal size
        of 80 characters regardless of the width of the window. (Default
        is None).
    """
    # Import some modules helpful for text formatting.
    from textwrap import TextWrapper
    from tabulate import tabulate
    from os import get_terminal_size

    # Try to get the actual number of columns in the terminal.
    if width is None:
        width = 80
        try:
            width = get_terminal_size().columns
        except Exception as e:
            logger.debug(f"Failed to get terminal size (using default "
                         f"{width}): {e}.")

        # Apply the maximum.
        if pretty_print_max_width is not None:
            assert isinstance(pretty_print_max_width, int)
            width = min(width, pretty_print_max_width)

    # Parameterize the text wrappers that format the ev text and the metadata.
    stmt_tr = TextWrapper(width=width)
    metadata_tr = TextWrapper(width=16)
    evidence_tr = TextWrapper(width=width - metadata_tr.width - 2)

    # Print the table.
    for i, s in enumerate(stmt_list[:stmt_limit]):

        # Print the Statement heading.
        stmt_str = f"[LIST INDEX: {i}] " + str(s)
        print(stmt_tr.fill(stmt_str))
        print("="*width)

        # Print the evidence
        for j, ev in enumerate(s.evidence[:ev_limit]):

            # Gather the metadata we want to display.
            metadata = [("EV INDEX", j), ("SOURCE", ev.source_api)]
            for id_type in ['PMID', 'PMCID', 'DOI']:
                if id_type in ev.text_refs:
                    metadata.append((id_type, ev.text_refs[id_type]))
                    break

            # Form the metadata string to fill out its allocated space.
            metadata_str = '\n'.join(line + ' '*(metadata_tr.width - len(line))
                                     for k, v in metadata
                                     for line in metadata_tr.wrap(f"{k}: {v}"))

            # Form the evidence string.
            if ev.text:
                text_str = evidence_tr.fill(ev.text)
            else:
                text_str = evidence_tr.fill("(No evidence text)")

            # Print the entire thing
            full_str = tabulate([[metadata_str, text_str]], tablefmt='plain')
            print(full_str)
            print('-'*width)
        print()


class UnresolvedUuidError(Exception):
    pass


class InputError(Exception):
    pass


