import pickle
from pysb import bng
from assemble_model import assemble_model
import pygraphviz

def get_braf_species(model):
    sp = [s for s in model.species if str(s).find('BRAF') != -1]
    sp = [s for s in sp if str(s).find('WT') == -1]
    return sp


def _get_node_name(species):
    name = ''
    monomer_codes = [m.monomer.name[0] for m in species.monomer_patterns]
    for mp in species.monomer_patterns:
        name += mp.monomer.name[0]
    return name


def draw_reaction_network(model, species):
    graph = pygraphviz.AGraph(directed=True, rankdir="LR")
    for sp in species:
        name = _get_node_name(sp)
        print(name)
        color = "#ccffcc"
        graph.add_node(name, label=name, shape='Mrecord', fillcolor=color,
                       style='filled', color='transparent', fontsize='12',
                       margin='0.06,0')
    species_set = set([model.species.index(s) for s in species])
    for i, r in enumerate(model.reactions_bidirectional):
        reactants = set(r['reactants'])
        products = set(r['products'])
        modifiers = reactants & products
        reactants = reactants - modifiers
        products = products - modifiers
        if not species_set.intersection(modifiers.union(reactants).union(products)):
            continue
        #reaction_node = 'r%d' % i
        #graph.add_node(reaction_node, label=reaction_node, shape='circle',
        #               fillcolor='lightgray', style='filled',
        #               color='transparent', fontsize='12',
        #               width=".3", height=".3", margin="0.06,0")
        attr_reversible = {'dir': 'both', 'arrowtail': 'empty'} if r['reversible'] else {}
        for s in reactants:
            if s not in species_set:
                continue
            for p in products:
                if p not in species_set:
                    continue
                source = _get_node_name(model.species[s])
                target = _get_node_name(model.species[p])
                graph.add_edge(source, target, **attr_reversible)
    return graph

if __name__ == '__main__':
    #model = assemble_model(4, False)
    #bng.generate_equations(model)
    with open('model_rxn_cache.pkl', 'rb') as fh:
        model = pickle.load(fh)
    braf_species = get_braf_species(model)
    graph = draw_reaction_network(model, braf_species)
    with open('braf_rxn_network.dot', 'wt') as fh:
        fh.write(graph.string())
    print(len(braf_species))
