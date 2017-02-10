import sys
from indra import reach, trips
from indra.assemblers import GraphAssembler

from run_ras_pathway import process_trips, process_reach

def draw_graph(stmts):
    graphpr = {'rankdir': 'TD'}
    nodepr = {'fontsize': 12, 'shape': 'plaintext', 'margin': '0,0', 'pad': 0}
    ga = GraphAssembler(stmts, graph_properties=graphpr, node_properties=nodepr)
    ga.make_model()
    ga.save_dot('p90rsk_correction.dot')
    ga.save_pdf('p90rsk_correction.pdf')

if __name__ == '__main__':
    orig_txt = [ln.strip() for ln in open('ras_pathway.txt', 'rt').readlines()]
    correct_txt = [ln.strip() for ln in open('correction.txt', 'rt').readlines()]

    for ln in correct_txt:
        if ln.startswith('<'):
            remove_line = ln[2:]
            orig_txt.remove(remove_line)
        elif ln.startswith('>'):
            add_line = ln[2:]
            orig_txt.append(add_line)

    txt = ' '.join(orig_txt)

    if len(sys.argv) < 2:
        print('Reader not specified')
        sys.exit()
    reader = sys.argv[1]
    if reader == 'reach':
        print('Using REACH')
        stmts = process_reach(txt, True)
    elif reader == 'trips':
        print('Using TRIPS')
        stmts = process_trips(txt, True)

    for st in stmts:
        print '%s\t%s' % (st, st.evidence[0].text)

    draw_graph(stmts)

