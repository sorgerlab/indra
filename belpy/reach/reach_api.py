import os
import json
import tempfile
from belpy.java_vm import autoclass, JavaException
from processor import ReachProcessor

_nxml_fries_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                   '../../../../groups/surdeanu/nxml2fries/nxml2fries')

def process_nxml(file_name):
    base = os.path.basename(file_name)
    file_id = os.path.splitext(base)[0]
    tmp_dir = tempfile.mkdtemp()
    try:
        paper_reader = autoclass('edu.arizona.sista.bionlp.ReadPaper')
        paper_reader.main([file_name, tmp_dir, _nxml_fries_path])
    except JavaException:
        print 'Could not process file %s.' % file_name
        return None
    
    json_file_name = os.path.join(tmp_dir, file_id + '.uaz.events.json')
    return process_json_file(json_file_name)

def process_json_file(file_name):
    try:
        with open(file_name, 'rt') as fh:
            json_str = fh.read()
            return process_json_str(json_str)
    except IOError:
        print 'Could not read file %s.' % file_name


def process_json_str(json_str):
    json_dict = json.loads(json_str)
    rp = ReachProcessor(json_dict)
    rp.get_phosphorylation()
    rp.get_complexes()
    return rp

if __name__ == '__main__':
    rp = process_json_file('PMC0000001.uaz.events.json')
