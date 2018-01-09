"""Provide an api used to run and get statements from the sparser reading tool.
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_version', 'process_text', 'process_xml', 'process_nxml_file',
           'process_nxml_str', 'make_sparser_nxml_from_text', 'run_sparser',
           'process_json_dict']

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
    assert sparser_path is not None, "Sparser path is not defined."
    with open(os.path.join(sparser_path, 'version.txt'), 'r') as f:
        version = f.read().strip()
    return version


def process_xml(xml_str):
    try:
        tree = ET.XML(xml_str, parser=UTB())
    except ET.ParseError as e:
        logger.error('Could not parse XML string')
        logger.error(e)
        return None
    sp = _process_elementtree(tree)
    return sp


def run_sparser(fname, output_fmt, outbuf=None):
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
    output_path = fname.split('.')[0] + '-semantics' + suffix
    for fpath in [sparser_exec_path, fname]:
        if not os.path.exists(fpath):
            raise Exception("'%s' is not a valid path." % fpath)

    out_bts = subprocess.check_output([sparser_exec_path, format_flag, fname])
    if outbuf is not None:
        outbuf.write(out_bts)
        outbuf.flush()
    assert os.path.exists(output_path), 'No output file created by sparser.'
    return output_path


def convert_sparser_output(output_fname, output_fmt='json'):
    if output_fmt not in ['json', 'xml']:
        logger.error("Unrecognized output format '%s'." % output_fmt)
        return None

    ret = None
    with open(output_fname, 'rt') as fh:
        if output_fmt == 'json':
            json_dict = json.load(fh)
            ret = process_json_dict(json_dict)
        else:
            xml_str = fh.read()
            ret = process_xml(xml_str)
    return ret


def process_nxml_file(fname, output_fmt='json', outbuf=None, cleanup=True):
    ret = None
    out_fname = None
    try:
        out_fname = run_sparser(fname, output_fmt, outbuf)
        ret = convert_sparser_output(out_fname, output_fmt)
    except Exception as e:
        logger.error("Sparser failed to run on %s." % fname)
        logger.exception(e)
    finally:
        if out_fname is not None and os.path.exists(out_fname) and cleanup:
            os.remove(out_fname)

    return ret


def make_sparser_nxml_from_text(text):
    text = _escape_xml(text)
    header = '<?xml version="1.0" encoding="UTF-8" ?>' + \
        '<OAI-PMH><article><body><sec id="s1"><p>'
    footer = '</p></sec></body></article></OAI-PMH>'
    nxml_str = header + text + footer
    return nxml_str


def process_text(text, output_fmt='json', outbuf=None, cleanup=True, key=''):
    nxml_str = make_sparser_nxml_from_text(text)
    return process_nxml_str(nxml_str, output_fmt, outbuf, cleanup, key)


def process_nxml_str(nxml_str, output_fmt='json', outbuf=None, cleanup=True,
                     key=''):
    tmp_fname = 'PMC%s_%d.nxml' % (key, mp.current_process().pid)
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
