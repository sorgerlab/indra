import urllib2
from functools32 import lru_cache
import xml.etree.ElementTree as et

hgnc_url = 'http://rest.genenames.org/fetch/'

@lru_cache(maxsize=1000)
def get_hgnc_name(hgnc_id):
    xml_tree = get_hgnc_entry(hgnc_id)
    if xml_tree is None:
        return None
    hgnc_name_tag =\
        xml_tree.find("result/doc/str[@name='symbol']")
    if hgnc_name_tag is None:
        return None
    return hgnc_name_tag.text.strip()

def get_hgnc_entry(hgnc_id):
    url = hgnc_url + 'hgnc_id/%s' % hgnc_id
    headers = {'Accept': '*/*'}
    req = urllib2.Request(url, headers=headers)
    try:
        res = urllib2.urlopen(req)
    except urllib2.HTTPError:
        return None
    xml_tree = et.parse(res)
    return xml_tree
