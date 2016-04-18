import urllib
import urllib2
import xml.etree.ElementTree as et
from indra.literature import pubmed_client

pmc_url = 'http://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi'

def get_ids(search_term, retmax=1000):
    return pubmed_client.get_ids(search_term, retmax=retmax, db='pmc')

def get_xml(pmc_id):
    if pmc_id.upper().startswith('PMC'):
        pmc_id = pmc_id[3:]

    params = {}
    params['verb'] = 'GetRecord'
    params['identifier'] = 'oai:pubmedcentral.nih.gov:%s' % pmc_id
    params['metadataPrefix'] = 'pmc'

    data = urllib.urlencode(params)
    req = urllib2.Request(pmc_url, data)
    try:
        res = urllib2.urlopen(req)
    except urllib2.HTTPError:
        print 'Couldn\'t download PMC%d' % pmc_id
    xml_str = res.read()
    xml_str = xml_str.decode('utf-8')

    err = check_xml_error(xml_str)
    if err is None:
        return xml_str
    else:
        print 'PMC client returned with error %s: %s' % (err[0], err[1])
        return None

def check_xml_error(xml_str):
    tree = et.fromstring(xml_str.encode('utf-8'))
    xmlns = "http://www.openarchives.org/OAI/2.0/"
    err_tag = tree.find('{%s}error' % xmlns)
    if err_tag is not None:
        err_code = err_tag.attrib['code']
        err_text = err_tag.text
        return (err_code, err_text)
    return None
