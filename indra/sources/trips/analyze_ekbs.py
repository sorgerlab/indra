import sys
import glob
import networkx
import xml.etree.ElementTree as ET


def build_event_graph(graph, tree, node):
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


def node_key(term):
    return term.attrib.get('id')


def get_text(tag):
    tt = tag.find('text')
    if tt is not None:
        return tt.text


def get_args(node):
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
    type = node.find('type')
    if type is not None:
        return type.text
    return None


def get_node_by_id(tree, node_id):
    for node_type in ['EVENT', 'TERM', 'CC']:
        node = tree.find("%s[@id='%s']" % (node_type, node_id))
        if node is not None:
            return node


def print_subtree(node, level):
    s = ''
    pre = '--' * level
    type = get_type(node)
    text = get_text(node)
    s += pre + type + ' | ' + text + '\n'
    args = get_args(node)
    for role, id in args.items():
        sub_node = get_node_by_id(et, id)
        if sub_node is not None:
            text = get_text(sub_node)
            s += pre + role + ' | ' + text + '\n'
            s += print_subtree(sub_node, level+1)
    return s


def type_match(a, b):
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
    ag = networkx.nx_agraph.to_agraph(graph)
    ag.draw(fname, prog='dot')


def build_patterns(fnames):
    patterns = []
    for fn in fnames:
        et = ET.parse(fn)
        res = et.findall('CC') + et.findall('EVENT')
        for event in res:
            G = networkx.DiGraph()
            build_event_graph(G, et, event)
            add_graph(patterns, G)
    patterns = sorted(patterns, key=lambda x: len(x), reverse=True)
    return patterns


if __name__ == '__main__':
    search_folder = sys.argv[1]
    fnames = glob.glob('%s/*.ekb' % search_folder)
    patterns = build_patterns(fnames)

    """
    cc_types = {}
    event_types = {}
    trees = {}
    for fn in fnames:
        et = ET.parse(fn)
        res = et.findall('CC')
        for event in res:
            type = get_type(event)
            if type is not None:
                try:
                    cc_types[type].append(event.attrib['id'])
                except KeyError:
                    cc_types[type] = [event.attrib['id']]
            s = print_subtree(event, 0)
            try:
                trees[s].append(event.attrib['id'])
            except KeyError:
                trees[s] = [event.attrib['id']]

        res = et.findall('EVENT')
        for event in res:
            type = get_type(event)
            if type is not None:
                try:
                    event_types[type].append(event.attrib['id'])
                except KeyError:
                    event_types[type] = [event.attrib['id']]
            s = print_subtree(event, 0)
            try:
                trees[s].append(event.attrib['id'])
            except KeyError:
                trees[s] = [event.attrib['id']]
    """
"""
class Pattern(object):
    def __init__(self, name, nodes, edges):
        self.name = name
        self.graph = networkx.DiGraph()
        for node_id, node_attrs in nodes.items():
            self.graph.add_node(node_id, **node_attrs)
        for source_id, out_edges in edges.items():
            for target_id, edge_attrs in out_edges.items():
                self.graph.add_edge(source_id, target_id, **edge_attrs)

    def draw(self, fname):
        ag = networkx.nx_agraph.to_agraph(self.graph)
        ag.draw(fname, prog='dot')

    def matches(self, other):
        if networkx.is_isomorphic(self.graph, other):
            return True
        return False

patterns = []

#####
nodes = {'factor': {'type': 'ONT::PROTEIN'},
         'cc': {'type': 'ONT::CAUSE'},
         'outcome': {'type': 'ONT::BIOLOGICAL-PROCESS'}}
edges = {'factor': {},
         'cc': {'factor': {'type': ':FACTOR'}, 'outcome': {'type': ':OUTCOME'}},
         'outcome': {}}
reg_protein_process = Pattern('reg_protein_process', nodes, edges)
patterns.append(reg_protein_process)
#####

def node_matches(n1, n2):
    return n1['type'] == n2['type']

def edge_matches(e1, e2):
    return e1['type'] == e2['type']

"""

