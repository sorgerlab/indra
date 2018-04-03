import logging
import xml.etree.ElementTree as ET
from collections import defaultdict, Counter
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('csxml')


def process_file(filename):
    logger.info("Parsing %s to XML" % filename)
    with open(filename, 'rb') as f:
        et = ET.parse(f, parser=UTB())
    return CsxmlProcessor(et)

class CsxmlProcessor(object):
    def __init__(self, tree):
        self._tree = tree
        self.matches = [(m.get('type'), m.get('name'), m.get('urn'),
                         m.get('msid'))
                        for m in tree.findall('./doc/sec/sent/match/entity')]
        self.match_ctr = Counter([t[0] for t in self.matches])

    # For each document
        # For each 
