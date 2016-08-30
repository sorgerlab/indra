import sys
import pickle
from indra import reach
from indra.assemblers import GraphAssembler

txt = open('extension.txt', 'rt').read()
rp = reach.process_text(txt, offline=True)
st = rp.statements
for s in st:
    print '%s\t%s' % (s, s.evidence[0].text)

with open('extension.pkl', 'wb') as fh:
    pickle.dump(st, fh)

graphpr = {'rankdir': 'TD'}
nodepr = {'fontsize': 12, 'shape': 'plaintext', 'margin': '0,0', 'pad': 0}
ga = GraphAssembler(st, graph_properties=graphpr, node_properties=nodepr)
ga.make_model()
ga.save_dot('jnk_extension.dot')
ga.save_pdf('jnk_extension.pdf')
