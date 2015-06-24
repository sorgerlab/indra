import urllib
import urllib2
from BeautifulSoup import BeautifulSoup

hgnc_url = 'http://www.genenames.org/cgi-bin/gene_symbol_report'

def get_hgnc_name(hgnc_id):
    html = get_hgnc_entry(hgnc_id)
    parsed_html = BeautifulSoup(html)
    hgnc_name = parsed_html.body.find('dd', attrs={'id':'app_symbol'}).text
    return hgnc_name

def get_hgnc_entry(hgnc_id):
    hgnc_args = urllib.urlencode({'hgnc_id': 'HGNC:%s' % hgnc_id})
    req = urllib2.Request(hgnc_url, hgnc_args)
    res = urllib2.urlopen(req)
    html = res.read()
    return html
