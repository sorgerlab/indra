import glob
import shutil
import json
import subprocess
import time
import sys
from indra.pysb_assembler import PysbAssembler
from indra.reach import reach_api
from indra.trips import trips_api
from reach_graph_assembler import ReachGraphAssembler

def process_text(path):
    fnames = glob.glob(path + '.txt')
    for fname in fnames:
        with open(fname, 'rt') as fh:
            print 'Processing %s' % fname
            txt = fh.read()
            rp = reach_api.process_text(txt)
            shutil.copyfile('reach_output.json', fname[:-4]+'.json')

def assemble_indra(path, model_name):
    fnames = glob.glob(path + '.json')
    pa = PysbAssembler()
    for fname in fnames:
        with open(fname, 'rt') as fh:
            print 'Processing %s' % fname
            json_str = fh.read()
            rp = reach_api.process_json_str(json_str, events_only=False)
            pa.add_statements(rp.statements)
    pa.make_model()
    pa.print_model(model_name + '.py')

def assemble_graph(path, model_name):
    file_names = glob.glob(path + '.json')
    ga = ReachGraphAssembler()
    for f in file_names:
        fh = open(f, 'rt')
        txt = fh.read()
        txt = txt.replace('frame-id', 'frame_id')
        txt = txt.replace('argument-label', 'argument_label')

        json_dict = json.loads(txt)
        print f, len(json_dict['events']['frames'])
        fh.close()
        
        ga.extend_graph(json_dict['events'], json_dict['entities'])
    
    with open(model_name + '.dot', 'wt') as fh:
        fh.write(ga.get_string())

    subprocess.call(('dot -T ps2 -o %s.ps %s.dot' % (model_name, model_name)).split(' '))
    subprocess.call(('dot -T png -o %s.png %s.dot' % (model_name, model_name)).split(' '))
    subprocess.call(('ps2pdf %s.ps' % model_name).split(' '))

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        folder = sys.argv[1]
        if len(sys.argv) >= 3:
            do_process_text = sys.argv[2]
        else:
            do_process_text = False
    else:
        folder = 'blog'
    ts = time.time()
    if do_process_text:
        process_text(path=folder+'/*')
    assemble_graph(folder+'/*', folder)
    te = time.time()
    print 'Time taken: %d' % (te-ts)
