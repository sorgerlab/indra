import rdflib
from indra.statements import Agent, Influence
from indra.assemblers import GraphAssembler

prefixes = """
    PREFIX causal: <http://worldmodelers.com/CauseEffect#>
    PREFIX cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
    """


def load_graph(fname):
    g = rdflib.Graph()
    with open(fname, 'rb') as fh:
        print('Started loading graph')
        g.parse(fh, format='json-ld')
        print('Finished loading graph')
    return g


if __name__ == '__main__':
    g = load_graph('cag.json-ld')
    query = prefixes + """
        SELECT ?rel
            ?causetext
            ?effecttext
            ?evtext
        WHERE {
            ?rel a causal:CausalAssertion .
            ?rel causal:has_cause ?cause .
            ?rel causal:has_effect ?effect .
            ?rel cco:has_text_value ?evtext .
            ?cause cco:has_text_value ?causetext .
            ?effect cco:has_text_value ?effecttext .
            }
        """
    res = g.query(query)
    rdict = {}
    for rel, causetext, effecttext, evtext in res:
        relid = rel.rsplit('#')[1]
        if relid not in rdict:
            rdict[relid] = [(causetext, effecttext)]
        else:
            rdict[relid].append((causetext, effecttext))
    stmts = []
    for relid, ces in rdict.items():
        cause = sorted(set([c[0] for c in ces]), key=lambda x: len(x))[0]
        effect = sorted(set([c[1] for c in ces]), key=lambda x: len(x))[0]
        print('%s -> %s' % (cause, effect))
        stmts.append(Influence(Agent(cause), Agent(effect)))
    ga = GraphAssembler(stmts)
    ga.make_model()
    ga.print_pdf('bbn_cag.pdf')

