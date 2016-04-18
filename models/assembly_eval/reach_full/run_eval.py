import os
import shutil
from indra.reach import reach_api
from indra.literature import pmc_client
from indra.preassembler.hierarchy_manager import entity_hierarchy as eh
from indra.preassembler.hierarchy_manager import modification_hierarchy as mh
from indra.preassembler import Preassembler, render_stmt_graph
from indra.pysb_assembler import PysbAssembler

def have_file(fname):
    return os.path.exists(fname)

def print_stmts(stmts, file_name):
    with open(file_name, 'wt') as fh:
        for s in stmts:
            agents = s.agent_list()
            db_refs = [('%s(%s)' % (a.name, a.db_refs)) 
                        for a in agents if a is not None]
            db_refs_str = ', '.join(db_refs)
            fh.write('%s\t%s\t%s\t%s\n' %
                     (s, db_refs_str, 'PMC'+s.evidence[0].pmid,
                      s.evidence[0].text))

if __name__ == '__main__':
    pmc_ids = ['PMC1234335', 'PMC3178447', 'PMC3690480',
               'PMC4345513', 'PMC534114']
    rerun = False

    pa = Preassembler(eh, mh)

    for pi in pmc_ids:
        print 'Reading %s...' % pi
        # If REACH already processed it then don't run it again
        if rerun or not have_file(pi + '.json'):
            if have_file(pi + '.txt'):
                txt = open(pi + '.txt').read()
                rp = reach_api.process_text(txt)
            elif have_file(pi + '.nxml'):
                rp = reach_api.process_nxml(pi + '.nxml')
            else:
                rp = reach_api.process_pmc(pi, save=True)
            shutil.move('reach_output.json', pi + '.json')
        else:
            rp = reach_api.process_json_file(pi + '.json')

        print '%s statements collected.' % len(rp.statements)
        pa.add_statements(rp.statements)

    print '%d statements collected in total.' % len(pa.stmts)
    duplicate_stmts = pa.combine_duplicates()
    print '%d statements after combining duplicates.' % len(duplicate_stmts)
    related_stmts = pa.combine_related()
    print '%d statements after combining related.' % len(related_stmts)

    # Print the statement graph
    graph = render_stmt_graph(related_stmts)
    graph.draw('reach_graph.pdf', prog='dot')
    # Print statement diagnostics
    print_stmts(pa.stmts, 'reach_statements.tsv')
    print_stmts(related_stmts, 'reach_related_statements.tsv')

    pya = PysbAssembler()
    pya.add_statements(related_stmts)
    model = pya.make_model()

    print 'PySB model has %d monomers and %d rules' %\
        (len(model.monomers), len(model.rules))
