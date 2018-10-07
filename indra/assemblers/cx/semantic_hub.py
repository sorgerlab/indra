import networkx
from collections import defaultdict

def find_hub(graph):
    pass


def get_aspect(cx, aspect_name):
    for entry in cx:
        if list(entry.keys())[0] == aspect_name:
            return entry[aspect_name]


def classify_nodes(graph, hub):
    node_stats = defaultdict(lambda: defaultdict(list))
    for u, v, data in graph.edges(data=True):
        if hub == u:
            h, o = u, v
            if data['i'] != 'Complex':
                node_stats[o]['up'].append(1)
            else:
                node_stats[o]['up'].append(0)
        elif hub == v:
            h, o = v, u
            if data['i'] != 'Complex':
                node_stats[o]['up'].append(-1)
            else:
                node_stats[o]['up'].append(0)
        else:
            continue
        node_stats[o]['interaction'].append(data['i'])

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


def get_node_by_name(graph, name):
    for id, attrs in graph.nodes(data=True):
        if attrs['n'] == name:
            return id


if __name__ == '__main__':
    import json
    with open('CDK13.cx', 'r') as fh:
        cx = json.load(fh)
    graph = cx_to_networkx(cx)

    node_classes = classify_nodes(graph, get_node_by_name(graph, 'CDK13'))
