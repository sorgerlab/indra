import sys
import pickle
from indra import trips
from indra import reach
from indra.assemblers import GraphAssembler


def process_reach(txt, reread):
    if reread:
        rp = reach.process_text(txt, offline=True)
        st = rp.statements
    else:
        rp = reach.process_json_file('reach_output.json')
        st = rp.statements
    for s in st:
        print('%s\t%s' % (s, s.evidence[0].text))
    return st

def process_trips(txt, reread):
    if reread:
        stmts = []
        sentences = txt.strip().split('\n')
        for sentence in sentences:
            print(sentence)
            tp = trips.process_text(txt)
            stmts += tp.statements
    else:
        tp = trips.process_xml(open('trips_output.xml', 'r').read())
        stmts = tp.statements
    for st in stmts:
        print('%s\t%s' % (st, st.evidence[0].text))
    return stmts

def draw_graph(stmts):
    graphpr = {'rankdir': 'TD'}
    nodepr = {'fontsize': 12, 'shape': 'plaintext', 'margin': '0,0', 'pad': 0}
    ga = GraphAssembler(stmts, graph_properties=graphpr,
                        node_properties=nodepr)
    ga.make_model()
    ga.save_dot('ras_pathway.dot')
    ga.save_pdf('ras_pathway.pdf')

if __name__ == '__main__':
    reread = True
    txt = open('ras_pathway.txt', 'rt').read()
    print('-----')
    print(txt)
    print('-----')
    if len(sys.argv) < 2:
        print('Reader not specified')
        sys.exit()
    reader = sys.argv[1]
    if reader == 'reach':
        print('Using REACH')
        stmts = process_reach(txt, reread)
    elif reader == 'trips':
        print('Using TRIPS')
        stmts = process_trips(txt, reread)

