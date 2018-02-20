"""Provides an API used to run and get Statements from the Sparser
reading system.
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['process_text', 'process_nxml_str', 'process_nxml_file',
           'process_sparser_output', 'process_json_dict', 'process_xml',
           'run_sparser', 'get_version', 'make_nxml_from_text']

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


def process_text(text, output_fmt='json', outbuf=None, cleanup=True, key=''):
    """Return processor with Statements extracted by reading text with Sparser.

    Parameters
    ----------
    text : str
        The text to be processed
    output_fmt: Optional[str]
        The output format to obtain from Sparser, with the two options being
        'json' and 'xml'. Default: 'json'
    outbuf : Optional[file]
        A file like object that the Sparser output is written to.
    cleanup : Optional[bool]
        If True, the temporary file created, which is used as an input
        file for Sparser, as well as the output file created by Sparser
        are removed. Default: True
    key : Optional[str]
        A key which is embedded into the name of the temporary file
        passed to Sparser for reading. Default is empty string.

    Returns
    -------
    SparserXMLProcessor or SparserJSONProcessor depending on what output
    format was chosen.
    """
    nxml_str = make_nxml_from_text(text)
    return process_nxml_str(nxml_str, output_fmt, outbuf, cleanup, key)


def process_nxml_str(nxml_str, output_fmt='json', outbuf=None, cleanup=True,
                     key=''):
    """Return processor with Statements extracted by reading an NXML string.

    Parameters
    ----------
    nxml_str : str
        The string value of the NXML-formatted paper to be read.
    output_fmt: Optional[str]
        The output format to obtain from Sparser, with the two options being
        'json' and 'xml'. Default: 'json'
    outbuf : Optional[file]
        A file like object that the Sparser output is written to.
    cleanup : Optional[bool]
        If True, the temporary file created in this function,
        which is used as an input file for Sparser, as well as the
        output file created by Sparser are removed. Default: True
    key : Optional[str]
        A key which is embedded into the name of the temporary file
        passed to Sparser for reading. Default is empty string.

    Returns
    -------
    SparserXMLProcessor or SparserJSONProcessor depending on what output
    format was chosen.
    """
    tmp_fname = 'PMC%s_%d.nxml' % (key, mp.current_process().pid)
    with open(tmp_fname, 'wb') as fh:
        fh.write(nxml_str.encode('utf-8'))
    try:
        sp = process_nxml_file(tmp_fname, output_fmt, outbuf, cleanup)
    finally:
        if cleanup and os.path.exists(tmp_fname):
            os.remove(tmp_fname)
    return sp


def process_nxml_file(fname, output_fmt='json', outbuf=None, cleanup=True):
    """Return processor with Statements extracted by reading an NXML file.

    Parameters
    ----------
    fname : str
        The path to the NXML file to be read.
    output_fmt: Optional[str]
        The output format to obtain from Sparser, with the two options being
        'json' and 'xml'. Default: 'json'
    outbuf : Optional[file]
        A file like object that the Sparser output is written to.
    cleanup : Optional[bool]
        If True, the output file created by Sparser is removed.
        Default: True

    Returns
    -------
    sp : SparserXMLProcessor or SparserJSONProcessor depending on what output
    format was chosen.
    """
    sp = None
    out_fname = None
    try:
        out_fname = run_sparser(fname, output_fmt, outbuf)
        sp = process_sparser_output(out_fname, output_fmt)
    except Exception as e:
        logger.error("Sparser failed to run on %s." % fname)
        logger.exception(e)
    finally:
        if out_fname is not None and os.path.exists(out_fname) and cleanup:
            os.remove(out_fname)

    return sp


def process_sparser_output(output_fname, output_fmt='json'):
    """Return a processor with Statements extracted from Sparser XML or JSON

    Parameters
    ----------
    output_fname : str
        The path to the Sparser output file to be processed. The file can
        either be JSON or XML output from Sparser, with the output_fmt
        parameter defining what format is assumed to be processed.
    output_fmt : Optional[str]
        The format of the Sparser output to be processed, can either be
        'json' or 'xml'. Default: 'json'

    Returns
    -------
    sp : SparserXMLProcessor or SparserJSONProcessor depending on what output
    format was chosen.
    """
    if output_fmt not in ['json', 'xml']:
        logger.error("Unrecognized output format '%s'." % output_fmt)
        return None

    sp = None
    with open(output_fname, 'rt') as fh:
        if output_fmt == 'json':
            json_dict = json.load(fh)
            sp = process_json_dict(json_dict)
        else:
            xml_str = fh.read()
            sp = process_xml(xml_str)
    return sp


def process_json_dict(json_dict):
    """Return processor with Statements extracted from a Sparser JSON.

    Parameters
    ----------
    json_dict : dict
        The JSON object obtained by reading content with Sparser, using the
        'json' output mode.

    Returns
    -------
    sp : SparserJSONProcessor
        A SparserJSONProcessor which has extracted Statements as its
        statements attribute.
    """
    sp = SparserJSONProcessor(json_dict)
    sp.get_statements()
    return sp


def process_xml(xml_str):
    """Return processor with Statements extracted from a Sparser XML.

    Parameters
    ----------
    xml_str : str
        The XML string obtained by reading content with Sparser, using the
        'xml' output mode.

    Returns
    -------
    sp : SparserXMLProcessor
        A SparserXMLProcessor which has extracted Statements as its
        statements attribute.
    """
    try:
        tree = ET.XML(xml_str, parser=UTB())
    except ET.ParseError as e:
        logger.error('Could not parse XML string')
        logger.error(e)
        return None
    sp = _process_elementtree(tree)
    return sp


def run_sparser(fname, output_fmt, outbuf=None):
    """Return the path to reading output after running Sparser reading.

    Parameters
    ----------
    fname : str
        The path to an input file to be processed. Due to the Spaser
        executable's assumptions, the file name needs to start with PMC
        and should be an NXML formatted file.
    output_fmt : Optional[str]
        The format in which Sparser should produce its output, can either be
        'json' or 'xml'.
    outbuf : Optional[file]
        A file like object that the Sparser output is written to.

    Returns
    -------
    output_path : str
        The path to the output file created by Sparser.
    """
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


def get_version():
    """Return the version of the Sparser executable on the path.

    Returns
    -------
    version : str
        The version of Sparser that is found on the Sparser path.
    """
    assert sparser_path is not None, "Sparser path is not defined."
    with open(os.path.join(sparser_path, 'version.txt'), 'r') as f:
        version = f.read().strip()
    return version


def make_nxml_from_text(text):
    """Return raw text wrapped in NXML structure.

    Parameters
    ----------
    text : str
        The raw text content to be wrapped in an NXML structure.

    Returns
    -------
    nxml_str : str
        The NXML string wrapping the raw text input.
    """
    text = _escape_xml(text)
    header = '<?xml version="1.0" encoding="UTF-8" ?>' + \
        '<OAI-PMH><article><body><sec id="s1"><p>'
    footer = '</p></sec></body></article></OAI-PMH>'
    nxml_str = header + text + footer
    return nxml_str


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
