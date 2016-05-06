import os
import sys
import pickle
import shutil
from indra import reach
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
    rerun = False

    # Download the papers if they are not available yet
    pmids = []
    for pmcid in pmc_ids:
        if not have_file(pmcid + '.nxml') and\
           not have_file(pmcid + '.txt'):
            txt, txt_format = get_full_text(pmcid)
            if txt_format == 'nxml':
                fname = pmcid + '.nxml'
            else:
                fname = pmcid + '.txt'
            with open(fname, 'wt') as fh:
                fh.write(txt.encode('utf-8'))
        pmids.append(id_lookup(pmcid)['pmid'])


    # Read each paper if it hasn't been read yet.
    # Otherwise use the existing json extractions.
    for pmcid, pmid in zip(pmc_ids, pmids):
        print 'Processing %s...' % pmcid
        # If REACH already processed it then don't run it again
        if rerun or not have_file(pmcid + '.json'):
            if have_file(pmcid + '.txt'):
                txt = open(pmcid + '.txt').read().decode('utf-8')
                rp = reach.process_text(txt, citation=pmid)
            elif have_file(pmcid + '.nxml'):
                rp = reach.process_nxml_file(pmcid + '.nxml', citation=pmid)
            shutil.move('reach_output.json', pmcid + '.json')
        else:
            rp = reach.process_json_file(pmcid + '.json', citation=pmid)

        # Instantiate the Preassembler
        pa = Preassembler(eh, mh)

        print '%d statements collected in total.' % len(pa.stmts)
        pa.add_statements(rp.statements)
        unique_stmts = pa.combine_duplicates()
        print '%d statements after combining duplicates.' % len(unique_stmts)
        related_stmts = pa.combine_related()
        print '%d statements after combining related.' % len(related_stmts)

        with open(pmcid + '.pkl', 'wb') as fh:
            pickle.dump(unique_stmts, fh)

        for st in sorted(unique_stmts, key=lambda x: len(x.evidence),
                                                     reverse=True):
            print len(st.evidence), st

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
