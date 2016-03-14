import os
import glob
import shutil
from indra.trips import trips_api
from indra.databases import pmc_client
from indra.preassembler.hierarchy_manager import entity_hierarchy as eh
from indra.preassembler.hierarchy_manager import modification_hierarchy as mh
from indra.preassembler import Preassembler
from indra.pysb_assembler import PysbAssembler

def have_file(fname):
    return os.path.exists(fname)

if __name__ == '__main__':
    fnames = glob.glob('*.ekb')
    
    pa = Preassembler(eh, mh)
    
    for fn in fnames:
        print '\n\n----------------------------'
        print 'Processing %s...' % fn
        xml_str = open(fn, 'rt').read()
        tp = trips_api.process_xml(xml_str)
        print 'Extracted events by type'
        print '------------------------'
        for k,v in tp.extracted_events.iteritems():
            print k, len(v)
        print '------------------------'
        print '%s statements collected.' % len(tp.statements)
        pa.add_statements(tp.statements)
        print '----------------------------\n\n'

    print '%d statements collected in total.' % len(pa.stmts)
    duplicate_stmts = pa.combine_duplicates()
    print '%d statements after combining duplicates.' % len(duplicate_stmts)
    related_stmts = pa.combine_related()
    print '%d statements after combining related.' % len(related_stmts)

    pya = PysbAssembler()
    pya.add_statements(related_stmts)
    model = pya.make_model()

    print 'PySB model has %d monomers and %d rules' %\
        (len(model.monomers), len(model.rules))
