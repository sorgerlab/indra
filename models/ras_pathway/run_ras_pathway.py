import sys
import pickle
from indra import reach
from indra.assemblers import GraphAssembler

if len(sys.argv) < 2:
    process_type = 'text'
else:
    process_type = sys.argv[1]

if process_type == 'text':
    txt = open('ras_pathway.txt', 'rt').read()
    rp = reach.process_text(txt, offline=True)
    st = rp.statements
elif process_type == 'json':
    rp = reach.process_json_file('reach_output.json')
    st = rp.statements
else:
    st = pickle.load(open('statements.pkl', 'rb'))
for s in st:
    print '%s\t%s' % (s, s.evidence[0].text)

graphpr = {'rankdir': 'TD'}
nodepr = {'fontsize': 12, 'shape': 'plaintext', 'margin': '0,0', 'pad': 0}
ga = GraphAssembler(st, graph_properties=graphpr, node_properties=nodepr)
ga.make_model()
ga.save_dot('ras_pathway.dot')
ga.save_pdf('ras_pathway.pdf')
