import json
import math
import random
import networkx
from collections import defaultdict

def find_hub(graph):
    pass


def get_aspect(cx, aspect_name):
    for entry in cx:
        if list(entry.keys())[0] == aspect_name:
            return entry[aspect_name]


def edge_type_to_class(edge_type):
    if 'Amount' in edge_type:
        return 'amount'
    if edge_type in ('Activation', 'Inhibition'):
        return 'activity'
    if edge_type == 'Complex':
        return 'complex'
    else:
        return 'modification'


def classify_nodes(graph, hub):
    node_stats = defaultdict(lambda: defaultdict(list))
    for u, v, data in graph.edges(data=True):
        # This means the node is downstream of the hub
        if hub == u:
            h, o = u, v
            if data['i'] != 'Complex':
                node_stats[o]['up'].append(-1)
            else:
                node_stats[o]['up'].append(0)
        # This means the node is upstream of the hub
        elif hub == v:
            h, o = v, u
            if data['i'] != 'Complex':
                node_stats[o]['up'].append(1)
            else:
                node_stats[o]['up'].append(0)
        else:
            continue
        node_stats[o]['interaction'].append(edge_type_to_class(data['i']))

    node_classes = {}
    for node_id, stats in node_stats.items():
        up = max(set(stats['up']), key=stats['up'].count)
        edge_type = max(set(stats['interaction']),
                        key=stats['interaction'].count)
        node_type = graph.nodes[node_id]['type']
        node_classes[node_id] = (up, edge_type, node_type)
    return node_classes


def get_attributes(aspect, id):
    attributes = {}
    for entry in aspect:
        if entry['po'] == id:
            attributes[entry['n']] = entry['v']
    return attributes


def cx_to_networkx(cx):
    graph = networkx.MultiDiGraph()
    for node_entry in get_aspect(cx, 'nodes'):
        id = node_entry['@id']
        attrs = get_attributes(get_aspect(cx, 'nodeAttributes'), id)
        attrs['n'] = node_entry['n']
        graph.add_node(id, **attrs)
    for edge_entry in get_aspect(cx, 'edges'):
        id = edge_entry['@id']
        attrs = get_attributes(get_aspect(cx, 'edgeAttributes'), id)
        attrs['i'] = edge_entry['i']
        graph.add_edge(edge_entry['s'], edge_entry['t'], key=id, **attrs)
    return graph


def get_quadrant_from_class(node_class):
    up, edge_type, _ = node_class
    if up == 0:
        return 0
    mappings = {(1, 'modification'): 1,
                (1, 'amount'): 2,
                (1, 'activity'): 3,
                (-1, 'activity'): 4,
                (-1, 'amount'): 5,
                (-1, 'modification'): 6}
    return mappings[(up, edge_type)]


def get_coordinates(node, node_class):
    quadrant_size = (2 * math.pi / 8.0)
    quadrant = get_quadrant_from_class(node_class)
    center_angle = (quadrant_size / 2.0) + quadrant_size * quadrant
    r = 100 + 200*random.random()
    alpha = center_angle - (quadrant_size / 2.0) * random.random() * quadrant_size
    x = r * math.cos(alpha)
    y = r * math.sin(alpha)
    return (x, y)


def get_layout_aspect(graph, hub, node_classes):
    aspect = []
    aspect.append({'@id': hub, 'x': 0.0, 'y': 0.0})
    for node, node_class in node_classes.items():
        if node == hub:
            continue
        x, y = get_coordinates(node, node_class)
        aspect.append({'@id': node, 'x': x, 'y': y})
    return aspect


def get_node_by_name(graph, name):
    for id, attrs in graph.nodes(data=True):
        if attrs['n'] == name:
            return id


if __name__ == '__main__':
    with open('CDK13.cx', 'r') as fh:
        cx = json.load(fh)
    graph = cx_to_networkx(cx)
    hub = 'CDK13'
    hub_node = get_node_by_name(graph, hub)
    node_classes = classify_nodes(graph, hub_node)
    layout_aspect = get_layout_aspect(graph, hub_node, node_classes)

