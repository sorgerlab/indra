import sys
import glob
import networkx
import xml.etree.ElementTree as ET
from indra.sources import trips


def node_key(term):
    """Return a key to be used for an element in the event EKB."""
    return term.attrib.get('id')


def get_text(tag):
    """Return the text associated with an elementin the event EKB."""
    tt = tag.find('text')
    if tt is not None:
        return tt.text


def get_args(node):
    """Return the arguments of a node in the event graph."""
    arg_roles = {}
    args = node.findall('arg') + \
        [node.find('arg1'), node.find('arg2'), node.find('arg3')]
    for arg in args:
        if arg is not None:
            id = arg.attrib.get('id')
            if id is not None:
                arg_roles[arg.attrib['role']] = (arg.attrib['id'], arg)
    # Now look at possible inevent links
    if node.find('features') is not None:
        inevents = node.findall('features/inevent')
        for inevent in inevents:
            if 'id' in inevent.attrib:
                arg_roles['inevent'] = (inevent.attrib['id'], inevent)

        ptms = node.findall('features/ptm') + node.findall('features/no-ptm')
        for ptm in ptms:
            if 'id' in inevent.attrib:
                arg_roles['ptm'] = (inevent.attrib['id'], ptm)

    # And also look for assoc-with links
    aw = node.find('assoc-with')
    if aw is not None:
        aw_id = aw.attrib['id']
        arg_roles['assoc-with'] = (aw_id, aw)
    return arg_roles


def get_type(node):
    """Return the type of an element in the event EKB."""
    type = node.find('type')
    if type is not None:
        return type.text
    return None


def type_match(a, b):
    """Return True of the types of a and b are compatible, False otherwise."""
    # If the types are the same, return True
    if a['type'] == b['type']:
        return True
    # Otherwise, look at some special cases
    eq_groups = [
        {'ONT::GENE-PROTEIN', 'ONT::GENE', 'ONT::PROTEIN'},
        {'ONT::PHARMACOLOGIC-SUBSTANCE', 'ONT::CHEMICAL'}
        ]
    for eq_group in eq_groups:
        if a['type'] in eq_group and b['type'] in eq_group:
            return True
    return False


def add_graph(patterns, G):
    """Add a graph to a set of unique patterns."""
    if not patterns:
        patterns.append([G])
        return
    for i, graphs in enumerate(patterns):
        if networkx.is_isomorphic(graphs[0], G, node_match=type_match,
                                  edge_match=type_match):
            patterns[i].append(G)
            return
    patterns.append([G])


def draw(graph, fname):
    """Draw a graph and save it into a file"""
    ag = networkx.nx_agraph.to_agraph(graph)
    ag.draw(fname, prog='dot')


def build_patterns(fnames):
    """Return a list of CC/EVENT graph patterns from a list of EKB files"""
    patterns = []
    for fn in fnames:
        et = ET.parse(fn)
        res = et.findall('CC') + et.findall('EVENT')
        for event in res:
            G = networkx.DiGraph()
            build_event_graph(G, et, event)
            add_graph(patterns, G)
    patterns = sorted(patterns, key=lambda x: len(x[0]), reverse=True)
    return patterns


def build_event_graph(graph, tree, node):
    """Return a DiGraph of a specific event structure, built recursively"""
    # If we have already added this node then let's return
    if node_key(node) in graph:
        return
    type = get_type(node)
    text = get_text(node)
    label = '%s (%s)' % (type, text)
    graph.add_node(node_key(node), type=type, label=label, text=text)
    args = get_args(node)
    for arg_role, (arg_id, arg_tag) in args.items():
        arg = get_node_by_id(tree, arg_id)
        if arg is None:
            arg = arg_tag
        build_event_graph(graph, tree, arg)
        graph.add_edge(node_key(node), node_key(arg), type=arg_role,
                       label=arg_role)


def get_extracted_events(fnames):
    """Get a full list of all extracted event IDs from a list of EKB files"""
    event_list = []
    for fn in fnames:
        tp = trips.process_xml_file(fn)
        ed = tp.extracted_events
        for k, v in ed.items():
            event_list += v
    return event_list


def check_event_coverage(patterns, event_list):
    """Calculate the ratio of patterns that were extracted."""
    proportions = []
    for pattern_list in patterns:
        proportion = 0
        for pattern in pattern_list:
            for node in pattern.nodes():
                if node in event_list:
                    proportion += 1.0 / len(pattern_list)
                    break
        proportions.append(proportion)
    return proportions


if __name__ == '__main__':
    search_folder = sys.argv[1]
    fnames = glob.glob('%s/*.ekb' % search_folder)
    patterns = build_patterns(fnames)
    extracted_events = get_extracted_events(fnames)
    proportions = check_event_coverage(patterns, extracted_events)
