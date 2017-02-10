import sys
import pickle
from indra import trips, reach
from indra.assemblers import GraphAssembler

def draw_graph(stmts):
    graphpr = {'rankdir': 'TD'}
    nodepr = {'fontsize': 12, 'shape': 'plaintext', 'margin': '0,0', 'pad': 0}
    ga = GraphAssembler(stmts, graph_properties=graphpr,
                        node_properties=nodepr)
    ga.make_model()
    ga.save_dot('jnk_extension.dot')
    ga.save_pdf('jnk_extension.pdf')

def process_reach(txt):
    rp = reach.process_text(txt, offline=True)
    return rp.statements

def process_trips(txt):
    tp = trips.process_text(txt)
    return tp.statements

if __name__ == '__main__':
    txt = open('extension.txt', 'rt').read()
    if len(sys.argv) < 2:
        print('Reader not specified')
        sys.exit()
    reader = sys.argv[1]
    if reader == 'reach':
        print('Using REACH')
        stmts = process_reach(txt)
    elif reader == 'trips':
        print('Using TRIPS')
        stmts = process_trips(txt)
    else:
        sys.exit()
    for s in stmts:
        print '%s\t%s' % (s, s.evidence[0].text)

    with open('extension.pkl', 'wb') as fh:
        pickle.dump(stmts, fh, protocol=2)

    draw_graph(stmts)
