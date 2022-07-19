from collections import Counter, defaultdict
import networkx
from indra.ontology.bio import bio_ontology


def plot_problem(problem):
    import matplotlib.pyplot as plt
    plt.ion()
    plt.figure()
    G = bio_ontology.subgraph(problem)
    pos = networkx.spring_layout(G)
    networkx.draw_networkx(G, pos, node_color='pink')
    edge_labels = networkx.get_edge_attributes(G, 'source')
    networkx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.show()


def print_cycle(cycle):
    for n1, n2 in zip(cycle, cycle[1:] + cycle[:1]):
        print('%s (%s)-[%s]-%s (%s)' % (
            bio_ontology.get_name(*bio_ontology.get_ns_id(n1)), n1,
            bio_ontology.edges[n1, n2]['type'],
            bio_ontology.get_name(*bio_ontology.get_ns_id(n2)), n2
            )
        )

if __name__ == '__main__':
    # First, find strongly connected components in the xref graph where
    # a given namespace appears more than once.
    bio_ontology.initialize()
    xrefs = [(e[0], e[1]) for e in bio_ontology.edges(data=True) if
             e[2]['type'] == 'xref']
    xrefg = bio_ontology.edge_subgraph(xrefs)
    comps = networkx.algorithms.strongly_connected_components(xrefg)

    problems = []
    for comp in comps:
        namespaces = [bio_ontology.get_ns(node) for node in comp]
        cnt = Counter(namespaces)
        if any(v > 1 for k, v in cnt.items()):
            problems.append(comp)

    print('Found %d problems in total' % len(problems))

    problems_by_ns = defaultdict(list)
    for problem in problems:
        nscnt = Counter([bio_ontology.get_ns(n) for n in problem])
        namespaces = [ns for ns, cnt in nscnt.items() if cnt > 1]
        for ns in namespaces:
            problems_by_ns[ns].append(problem)

    for ns, problems_ns in problems_by_ns.items():
        print(ns, len(problems_ns))


    # Next, find cycles in the isa/partof subgraph meaning circular
    # hierarchical relationships.
    hierarchy = [(e[0], e[1]) for e in bio_ontology.edges(data=True) if
                 e[2]['type'] in {'isa', 'partof'}]
    hierarchyg = bio_ontology.edge_subgraph(hierarchy)
    cycles = networkx.simple_cycles(hierarchyg)
    for cycle in cycles:
        print('---')
        print_cycle(cycle)
