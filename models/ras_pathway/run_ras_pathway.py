import sys
import time
import pickle
from indra.sources import reach, trips
from indra.assemblers import GraphAssembler
import indra.tools.assemble_corpus as ac


def process_reach(txt):
    print('Using REACH')
    ts = time.time()
    rp = reach.process_text(txt, offline=False)
    for s in rp.statements:
        print('%s\t%s' % (s, s.evidence[0].text))
    te = time.time()
    print('Time taken: %.2fs' % (te-ts))
    return rp.statements

def process_trips(txt, reread=True):
    print('Using TRIPS')
    ts = time.time()
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
    te = time.time()
    print('Time taken: %.2fs' % (te-ts))
    for st in stmts:
        print('%s\t%s' % (st, st.evidence[0].text))
    return stmts

def assemble_ras_pathway(fname, reader):
    # Make original pathway map
    with open(fname, 'rt') as fh:
        txt = fh.read()
    if reader == 'reach':
        stmts = process_reach(txt)
    elif reader == 'trips':
        stmts = process_trips(txt, reread=True)
    ac.dump_statements(stmts, 'ras_pathway.pkl')
    draw_graph(stmts, 'ras_pathway')

def assemble_correction(fname_orig, fname, reader):
    # Read correction
    with open(fname_orig, 'rt') as fh:
        orig_txt = [ln.strip() for ln in fh.readlines()]
    with open(fname, 'rt') as fh:
        correct_txt = [ln.strip() for ln in fh.readlines()]
    for ln in correct_txt:
        if ln.startswith('<'):
            remove_line = ln[2:]
            orig_txt.remove(remove_line)
        elif ln.startswith('>'):
            add_line = ln[2:]
            orig_txt.append(add_line)
    txt = '\n'.join(orig_txt)
    if reader == 'reach':
        stmts = process_reach(txt)
    elif reader == 'trips':
        stmts = process_trips(txt, reread=True)
    ac.dump_statements(stmts, 'ras_pathway_correction.pkl')
    draw_graph(stmts, 'ras_pathway_correction')

def assemble_extension(fname_orig, fname, reader):
    with open(fname_orig, 'rt') as fh:
        orig_txt = fh.read()
    with open(fname, 'rt') as fh:
        extension_txt = fh.read()
    txt = '\n'.join([orig_txt, extension_txt])
    if reader == 'reach':
        stmts = process_reach(txt)
    elif reader == 'trips':
        stmts = process_trips(txt, reread=True)
    ac.dump_statements(stmts, 'ras_pathway_extension.pkl')
    draw_graph(stmts, 'ras_pathway_extension')

def draw_graph(stmts, fname):
    graphpr = {'rankdir': 'TD'}
    nodepr = {'fontsize': 12, 'shape': 'plaintext', 'margin': '0,0', 'pad': 0}
    ga = GraphAssembler(stmts, graph_properties=graphpr,
                        node_properties=nodepr)
    ga.make_model()
    ga.save_dot('%s.dot' % fname)
    ga.save_pdf('%s.pdf' % fname)

if __name__ == '__main__':
    reread = True
    if len(sys.argv) < 2:
        print('Reader not specified')
        sys.exit()
    reader = sys.argv[1]
    if reader not in ['trips', 'reach']:
        print('Invaid reader')
        sys.exit()

    assemble_ras_pathway('ras_pathway.txt', reader)
    assemble_correction('ras_pathway.txt', 'correction.txt', reader)
    assemble_extension('ras_pathway.txt', 'extension.txt', reader)
