import rdflib
import collections
from indra.statements import Agent, Influence
from indra.sources.bbn import processor

def process_json_file(fname):
    return processor.BBNProcessor(fname)

if __name__ == '__main__':
    process_json_file('/Users/daniel/Downloads/bbn-m6-cag.v0.1/cag.json-ld')
