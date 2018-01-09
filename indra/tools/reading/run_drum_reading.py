import sys
import json
from indra.sources.trips.drum_reader import DrumReader
from indra.sources.trips import process_xml

def read_content(content):
    sentences = []
    for k, v in content.items():
        sentences += v
    dr = DrumReader(to_read=sentences)
    try:
        dr.start()
    except SystemExit:
        pass
    statements = []
    for extraction in dr.extractions:
        statements += process_xml(extraction).statements
    return statements

if __name__ == '__main__':
    file_name = sys.argv[1]
    with open(file_name, 'rt') as fh:
        content = json.load(fh)
    statements = read_content(content)
    print(statements)
