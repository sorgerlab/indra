import sys
from indra import reach
from indra.assemblers import GraphAssembler

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

rp = reach.process_text(txt, offline=True)
st = rp.statements
for s in st:
    print '%s\t%s' % (s, s.evidence[0].text)

graphpr = {'rankdir': 'TD'}
nodepr = {'fontsize': 12, 'shape': 'plaintext', 'margin': '0,0', 'pad': 0}
ga = GraphAssembler(st, graph_properties=graphpr, node_properties=nodepr)
ga.make_model()
ga.save_dot('rps6ka_correction.dot')
ga.save_pdf('rps6k1_correction.pdf')
