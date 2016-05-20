import re
import sys
import getopt
import urllib
import urllib2
import xml.dom.minidom
import logging

logger = logging.getLogger('trips')

trips_url = 'http://trips.ihmc.us/parser/cgi/drum'


def send_query(text, query_args=None):
    """Send a query to the TRIPS web service.

    Parameters
    ----------
    text : str
        The text to be processed.
    query_args : Optional[dict]
        A dictionary of arguments to be passed with the query.

    Returns
    -------
    html : str
        The HTML result returned by the web service.
    """
    if query_args is None:
        qa = {}
    qa['input'] = text
    data = urllib.urlencode(qa)
    req = urllib2.Request(trips_url, data)
    res = urllib2.urlopen(req, timeout=3600)
    html = res.read()
    return html


def get_xml(html):
    """Extract the EKB XML from the HTML output of the TRIPS web service.

    Parameters
    ----------
    html : str
        The HTML output from the TRIPS web service.

    Returns
    -------
    The extraction knowledge base (EKB) XML that contains the event and term
    extractions.
    """
    ekb = re.findall(r'<ekb.*?>(.*?)</ekb>', html, re.MULTILINE | re.DOTALL)
    if ekb:
        events_terms = ekb[0]
    else:
        events_terms = ''
    header = '<?xml version="1.0" encoding="utf-8" standalone="yes"?><ekb>'
    footer = '</ekb>'
    return header + events_terms + footer


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
        fh.write(xml_str_pretty.encode('utf8'))
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
            print 'String mode: python trips_client.py --string "RAS binds GTP" --output text.xml'
            print 'File mode: python trips_client.py --file test.txt --output text.xml'
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
            print 'Could not open %s.' % infile_name
            exit()
        text = fh.read()
        fh.close()
        print 'Parsing contents of %s...' % infile_name
    else:
        print 'Parsing string: %s' % text

    html = send_query(text)
    xml = get_xml(html)
    save_xml(xml, outfile_name)
