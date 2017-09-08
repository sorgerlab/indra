from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import json
import logging
import subprocess
import xml.etree.ElementTree as ET
from indra.util import UnicodeXMLTreeBuilder as UTB
from .processor import SparserXMLProcessor, SparserJSONProcessor

logger = logging.getLogger('sparser')

sparser_path_var = 'SPARSERPATH'
sparser_path = os.environ.get(sparser_path_var)

def process_xml(xml_str):
    try:
        tree = ET.XML(xml_str, parser=UTB())
    except ET.ParseError as e:
        logger.error('Could not parse XML string')
        logger.error(e)
        return None
    sp = _process_elementtree(tree)
    return sp

def process_nxml_file(fname, output_format='json'):
    if not sparser_path or not os.path.exists(sparser_path):
        logger.error('Sparser executable not set in %s' % sparser_path_var)
        return None
    if output_format == 'xml':
        format_flag = '-x'
        suffix = '.xml'
    elif output_format == 'json':
        format_flag = '-j'
        suffix = '.json'
    else:
        logger.error('Unknown output format: %s' % output_format)
        return None
    sparser_exec_path = os.path.join(sparser_path, 'save-semantics.sh')
    subprocess.call([sparser_exec_path, format_flag, fname])

    output_fname = fname.split('.')[0] + '-semantics' + suffix
    with open(output_fname, 'rt') as fh:
        if output_format == 'json':
            json_dict = json.load(fh)
            return process_json_dict(json_dict)
        else:
            xml_str = fh.read()
            return process_xml(xml_str)

def process_text(text, output_format='json'):
    header = '<?xml version="1.0" encoding="UTF-8" ?>' + \
        '<OAI-PMH><article><body><sec id="s1"><p>'
    footer = '</p></sec></body></article></OAI-PMH>'
    nxml_str = header + text + footer
    return process_nxml_str(nxml_str, output_format)


def process_nxml_str(nxml_str, output_format='json'):
    tmp_fname = 'PMC12345.nxml'
    with open(tmp_fname, 'wb') as fh:
        fh.write(nxml_str.encode('utf-8'))
    return process_nxml_file(tmp_fname, output_format)


def process_json_dict(json_dict):
    sp = SparserJSONProcessor(json_dict)
    sp.get_statements()
    return sp

def _process_elementtree(tree):
    sp = SparserXMLProcessor(tree)
    sp.get_modifications()
    sp.get_activations()
    return sp
