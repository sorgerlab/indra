from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import sys
import getopt
import xml.dom.minidom
import logging
import requests

logger = logging.getLogger('trips')

base_url = 'http://trips.ihmc.us/parser/cgi/'


def send_query(text, service_endpoint='drum', query_args=None):
    """Send a query to the TRIPS web service.

    Parameters
    ----------
    text : str
        The text to be processed.
    service_endpoint : Optional[str]
        Selects the TRIPS/DRUM web service endpoint to use. Is a choice between
        "drum" (default), "drum-dev", a nightly build, and "cwms" for use with
        more general knowledge extraction.
    query_args : Optional[dict]
        A dictionary of arguments to be passed with the query.

    Returns
    -------
    html : str
        The HTML result returned by the web service.
    """
    if service_endpoint in ['drum', 'drum-dev', 'cwms', 'cwmsreader']:
        url = base_url + service_endpoint
    else:
        logger.error('Invalid service endpoint: %s' % service_endpoint)
        return ''
    if query_args is None:
        query_args = {}
    query_args.update({'input': text})
    res = requests.get(url, query_args, timeout=3600)
    if not res.status_code == 200:
        logger.error('Problem with TRIPS query: status code %s' %
                     res.status_code)
        return ''
    # Gets unicode content
    return res.text


def get_xml(html, content_tag='ekb', fail_if_empty=False):
    """Extract the content XML from the HTML output of the TRIPS web service.

    Parameters
    ----------
    html : str
        The HTML output from the TRIPS web service.
    content_tag : str
        The xml tag used to label the content. Default is 'ekb'.
    fail_if_empty : bool
        If True, and if the xml content found is an empty string, raise an
        exception. Default is False.

    Returns
    -------
    The extraction knowledge base (e.g. EKB) XML that contains the event and
    term extractions.
    """
    cont = re.findall(r'<%(tag)s(.*?)>(.*?)</%(tag)s>' % {'tag': content_tag},
                      html, re.MULTILINE | re.DOTALL)
    if cont:
        events_terms = ''.join([l.strip() for l in cont[0][1].splitlines()])
        if 'xmlns' in cont[0][0]:
            meta = ' '.join([l.strip() for l in cont[0][0].splitlines()])
        else:
            meta = ''
    else:
        events_terms = ''
        meta = ''

    if fail_if_empty:
        assert events_terms != '',\
            "Got empty string for events content from html:\n%s" % html

    header = ('<?xml version="1.0" encoding="utf-8" standalone="yes"?><%s%s>'
              % (content_tag, meta))
    footer = '</%s>' % content_tag
    return header + events_terms.replace('\n', '') + footer


def save_xml(xml_str, file_name, pretty=True):
    """Save the TRIPS EKB XML in a file.

    Parameters
    ----------
    xml_str : str
        The TRIPS EKB XML string to be saved.
    file_name : str
        The name of the file to save the result in.
    pretty : Optional[bool]
        If True, the XML is pretty printed.
    """
    try:
        fh = open(file_name, 'wt')
    except IOError:
        logger.error('Could not open %s for writing.' % file_name)
        return
    if pretty:
        xmld = xml.dom.minidom.parseString(xml_str)
        xml_str_pretty = xmld.toprettyxml()
        fh.write(xml_str_pretty)
    else:
        fh.write(xml_str)
    fh.close()


if __name__ == '__main__':
    filemode = False
    text = 'Active BRAF phosphorylates MEK1 at Ser222.'
    outfile_name = 'braf_test.xml'

    opts, extraparams = getopt.getopt(sys.argv[1:], 's:f:o:h',
                                      ['string=', 'file=', 'output=', 'help'])
    for o, p in opts:
        if o in ['-h', '--help']:
            print('String mode: python trips_client.py '
                  '--string "RAS binds GTP" --output text.xml')
            print('File mode: python trips_client.py --file test.txt '
                  '--output text.xml')
            sys.exit()
        elif o in ['-s', '--string']:
            text = p
        elif o in ['-f', '--file']:
            filemode = True
            infile_name = p
        elif o in ['-o', '--output']:
            outfile_name = p

    if filemode:
        try:
            fh = open(infile_name, 'rt')
        except IOError:
            print('Could not open %s.' % infile_name)
            exit()
        text = fh.read()
        fh.close()
        print('Parsing contents of %s...' % infile_name)
    else:
        print('Parsing string: %s' % text)

    html = send_query(text)
    xml = get_xml(html)
    save_xml(xml, outfile_name)
