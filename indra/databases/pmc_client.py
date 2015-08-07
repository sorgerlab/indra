import urllib
import urllib2

def get_xml(pmc_id):
    pmc_url = 'http://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi'
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
    return xml_str
