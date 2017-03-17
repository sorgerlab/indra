import sys
import time
import pickle
from indra import trips
from indra import reach
from indra.assemblers import GraphAssembler
import indra.tools.assemble_corpus as ac


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
            tp = trips.process_text(sentence)
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
        ts = time.time()
        stmts = process_reach(txt, reread)
        te = time.time()
    elif reader == 'trips':
        print('Using TRIPS')
        ts = time.time()
        stmts = process_trips(txt, reread)
        te = time.time()
    else:
        print('Invalid reader')
        sys.exit()
    print('Time taken: %.2fs' % (te-ts))
    ac.dump_statements(stmts, 'statements.pkl')
