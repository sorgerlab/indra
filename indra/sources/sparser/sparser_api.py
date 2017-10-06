from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import json
import logging
import subprocess
import xml.etree.ElementTree as ET
import multiprocessing as mp
from indra.util import UnicodeXMLTreeBuilder as UTB
from .processor import SparserXMLProcessor, SparserJSONProcessor

logger = logging.getLogger('sparser')

sparser_path_var = 'SPARSERPATH'
sparser_path = os.environ.get(sparser_path_var)


def get_version():
    return '?'


def process_xml(xml_str):
    try:
        tree = ET.XML(xml_str, parser=UTB())
    except ET.ParseError as e:
        logger.error('Could not parse XML string')
        logger.error(e)
        return None
    sp = _process_elementtree(tree)
    return sp


def process_nxml_file(fname, output_fmt='json', outbuf=None, cleanup=True):
    if not sparser_path or not os.path.exists(sparser_path):
        logger.error('Sparser executable not set in %s' % sparser_path_var)
        return None
    if output_fmt == 'xml':
        format_flag = '-x'
        suffix = '.xml'
    elif output_fmt == 'json':
        format_flag = '-j'
        suffix = '.json'
    else:
        logger.error('Unknown output format: %s' % output_fmt)
        return None
    sparser_exec_path = os.path.join(sparser_path, 'save-semantics.sh')
    output_fname = fname.split('.')[0] + '-semantics' + suffix

    for fpath in [sparser_exec_path, fname]:
        if not os.path.exists(fpath):
            raise Exception("'%s' is not a valid path." % fpath)

    ret = None
    try:
        subprocess.call(
            [sparser_exec_path, format_flag, fname],
            stdout=outbuf,
            stdin=subprocess.PIPE
            )
        with open(output_fname, 'rt') as fh:
            if output_fmt == 'json':
                json_dict = json.load(fh)
                ret = process_json_dict(json_dict)
            else:
                xml_str = fh.read()
                ret = process_xml(xml_str)
    except Exception as e:
        logger.error("Sparser failed to run on %s." % fname)
        logger.exception(e)
    finally:
        if os.path.exists(output_fname) and cleanup:
            os.remove(output_fname)
    return ret


def process_text(text, output_fmt='json', outbuf=None, cleanup=True):
    text = _escape_xml(text)
    header = '<?xml version="1.0" encoding="UTF-8" ?>' + \
        '<OAI-PMH><article><body><sec id="s1"><p>'
    footer = '</p></sec></body></article></OAI-PMH>'
    nxml_str = header + text + footer
    return process_nxml_str(nxml_str, output_fmt, outbuf, cleanup)


def process_nxml_str(nxml_str, output_fmt='json', outbuf=None, cleanup=True):
    tmp_fname = 'PMC%d.nxml' % mp.current_process().pid
    with open(tmp_fname, 'wb') as fh:
        fh.write(nxml_str.encode('utf-8'))
    try:
        ret = process_nxml_file(tmp_fname, output_fmt, outbuf, cleanup)
    finally:
        if cleanup and os.path.exists(tmp_fname):
            os.remove(tmp_fname)
    return ret


def process_json_dict(json_dict):
    sp = SparserJSONProcessor(json_dict)
    sp.get_statements()
    return sp


def _process_elementtree(tree):
    sp = SparserXMLProcessor(tree)
    sp.get_modifications()
    sp.get_activations()
    return sp


def _escape_xml(text):
    esc_map = {'"': '&quot;', '&': '&amp;', '\'': '&apos;',
               '<': '&lt;', '>': '&gt;'}
    for orig, new in esc_map.items():
        text = text.replace(orig, new)
    return text
