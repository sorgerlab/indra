import os
import sys
import shutil
from indra import trips
from indra.literature import pmc_client, get_full_text, id_lookup
from indra.preassembler.hierarchy_manager import entity_hierarchy as eh
from indra.preassembler.hierarchy_manager import modification_hierarchy as mh
from indra.preassembler import Preassembler, render_stmt_graph
from indra.assemblers import PysbAssembler
from indra.assemblers import IndexCardAssembler

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

    pmids = [id_lookup(pmcid)['pmid'] for pmcid in pmc_ids]

    # Use the existing EKB extractions.
    for pmcid, pmid in zip(pmc_ids, pmids):
        print 'Processing %s...' % pmcid
        tp = trips.process_xml(open(pmcid + '-20160503T1152.ekb').read())

        # Instantiate the Preassembler
        pa = Preassembler(eh, mh)

        print '%s statements collected.' % len(tp.statements)
        pa.add_statements(tp.statements)
        print '%d statements collected in total.' % len(pa.stmts)
        duplicate_stmts = pa.combine_duplicates()
        print '%d statements after combining duplicates.' % len(duplicate_stmts)
        related_stmts = pa.combine_related()
        print '%d statements after combining related.' % len(related_stmts)

        # Assemble IndexCards
        ia = IndexCardAssembler(related_stmts)
        ia.make_model()
        ia.save_model(pmcid + '_card.json')
        # Print the statement graph
        graph = render_stmt_graph(related_stmts)
        graph.draw(pmcid + '_graph.pdf', prog='dot')
        # Print statement diagnostics
        print_stmts(pa.stmts, pmcid + '_statements.tsv')
        print_stmts(related_stmts, pmcid + '_related_statements.tsv')

        pya = PysbAssembler()
        pya.add_statements(related_stmts)
        model = pya.make_model()

        print 'PySB model has %d monomers and %d rules' %\
            (len(model.monomers), len(model.rules))
