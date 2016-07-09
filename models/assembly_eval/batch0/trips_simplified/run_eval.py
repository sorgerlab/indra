import os
import glob
import shutil
from indra import trips
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler import Preassembler, render_stmt_graph
from indra.assemblers import PysbAssembler

def have_file(fname):
    return os.path.exists(fname)

def print_stmts(stmts, file_name):
    with open(file_name, 'wt') as fh:
        for s in stmts:
            agents = s.agent_list()
            db_refs = [('%s(%s)' % (a.name, a.db_refs))
                        for a in agents if a is not None]
            db_refs_str = ', '.join(db_refs)
            fh.write('%s\t%s\t%s\n' %
                     (s, db_refs_str, s.evidence[0].text))

if __name__ == '__main__':
    fnames = glob.glob('*.txt')

    pa = Preassembler(hierarchies)

    for fn in fnames:
        print '\n\n----------------------------'
        print 'Processing %s...' % fn
        txt = open(fn, 'rt').read()
        tp = trips.process_text(txt)
        pa.add_statements(tp.statements)
        print '----------------------------\n\n'

    print '%d statements collected in total.' % len(pa.stmts)
    duplicate_stmts = pa.combine_duplicates()
    print '%d statements after combining duplicates.' % len(duplicate_stmts)
    related_stmts = pa.combine_related()
    print '%d statements after combining related.' % len(related_stmts)

    # Print the statement graph
    graph = render_stmt_graph(related_stmts)
    graph.draw('trips_graph.pdf', prog='dot')
    # Print statement diagnostics
    print_stmts(pa.stmts, 'trips_statements.tsv')
    print_stmts(related_stmts, 'trips_related_statements.tsv')

    pya = PysbAssembler()
    pya.add_statements(related_stmts)
    model = pya.make_model()

    print 'PySB model has %d monomers and %d rules' %\
        (len(model.monomers), len(model.rules))
