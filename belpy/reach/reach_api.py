import json

from reach_processor import ReachProcessor


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
    return rp

if __name__ == '__main__':
    rp = process_json_file('PMC0000001.uaz.events.json')
