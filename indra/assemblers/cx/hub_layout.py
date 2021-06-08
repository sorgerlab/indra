"""This module allows adding a semantic hub layout to NDEx CX networkx. This
is useful when a network is centered around a single hub node. The
layout generated here allocates different classes of nodes into segments
around the hub and then gives them random coordinates within that segment."""

import json
import math
import random
import networkx
from collections import defaultdict


def get_aspect(cx, aspect_name):
    """Return an aspect given the name of the aspect"""
    if isinstance(cx, dict):
        return cx.get(aspect_name)
    for entry in cx:
        if list(entry.keys())[0] == aspect_name:
            return entry[aspect_name]


def edge_type_to_class(edge_type):
    """Return the edge class for layout purposes based on the edge type"""
    edge_type = edge_type.lower()
    if 'amount' in edge_type:
        return 'amount'
    if edge_type in ('activation', 'inhibition'):
        return 'activity'
    if edge_type == 'complex':
        return 'complex'
    else:
        return 'modification'


def classify_nodes(graph, hub: int):
    """Classify each node based on its type and relationship to the hub."""
    node_stats = defaultdict(lambda: defaultdict(list))
    for u, v, data in graph.edges(data=True):
        # This means the node is downstream of the hub
        if hub == u:
            h, o = u, v
            if data['i'] != 'complex':
                node_stats[o]['up'].append(-1)
            else:
                node_stats[o]['up'].append(0)
        # This means the node is upstream of the hub
        elif hub == v:
            h, o = v, u
            if data['i'] != 'complex':
                node_stats[o]['up'].append(1)
            else:
                node_stats[o]['up'].append(0)
        else:
            continue
        node_stats[o]['interaction'].append(edge_type_to_class(data['i']))

    node_classes = {}
    for node_id, stats in node_stats.items():
        up = max(set(stats['up']), key=stats['up'].count)
        # Special case: if up is not 0 then we should exclude complexes
        # from the edge_type states so that we don't end up with
        # (-1, complex, ...) or (1, complex, ...) as the node class
        interactions = [i for i in stats['interaction'] if
                        not (up != 0 and i == 'complex')]
        edge_type = max(set(interactions), key=interactions.count)
        node_type = graph.nodes[node_id]['type']
        node_classes[node_id] = (up, edge_type, node_type)
    return node_classes


def get_attributes(aspect, id):
    """Return the attributes pointing to a given ID in a given aspect."""
    attributes = {}
    for entry in aspect:
        if entry['po'] == id:
            attributes[entry['n']] = entry['v']
    return attributes


def cx_to_networkx(cx):
    """Return a MultiDiGraph representation of a CX network."""
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
    """Return the ID of the segment of the plane corresponding to a class."""
    up, edge_type, _ = node_class
    if up == 0:
        return 0 if random.random() < 0.5 else 7
    mappings = {(-1, 'modification'): 1,
                (-1, 'amount'): 2,
                (-1, 'activity'): 3,
                (1, 'activity'): 4,
                (1, 'amount'): 5,
                (1, 'modification'): 6}
    return mappings[(up, edge_type)]


def get_coordinates(node_class):
    """Generate coordinates for a node in a given class."""
    quadrant_size = (2 * math.pi / 8.0)
    quadrant = get_quadrant_from_class(node_class)
    begin_angle = quadrant_size * quadrant
    r = 200 + 800*random.random()
    alpha = begin_angle + random.random() * quadrant_size
    x = r * math.cos(alpha)
    y = r * math.sin(alpha)
    return x, y


def get_layout_aspect(hub, node_classes):
    """Get the full layout aspect with coordinates for each node."""
    aspect = [{'node': hub, 'x': 0.0, 'y': 0.0}]
    for node, node_class in node_classes.items():
        if node == hub:
            continue
        x, y = get_coordinates(node_class)
        aspect.append({'node': node, 'x': x, 'y': y})
    return aspect


def get_node_by_name(graph, name):
    """Return a node ID given its name."""
    for id, attrs in graph.nodes(data=True):
        if attrs['n'] == name:
            return id


def add_semantic_hub_layout(cx, hub: str):
    """Attach a layout aspect to a CX network given a hub node."""
    graph = cx_to_networkx(cx)
    hub_node = get_node_by_name(graph, hub)
    node_classes = classify_nodes(graph, hub_node)
    layout_aspect = get_layout_aspect(hub_node, node_classes)
    cx['cartesianLayout'] = layout_aspect


if __name__ == '__main__':
    with open('CDK13.cx', 'r') as fh:
        cx = json.load(fh)
    add_semantic_hub_layout(cx, 'CDK13')
