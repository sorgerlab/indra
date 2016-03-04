import os
import csv
from functools32 import lru_cache
import urllib2
import xml.etree.ElementTree as et

hgnc_url = 'http://rest.genenames.org/fetch/'
hgnc_file = os.path.dirname(os.path.abspath(__file__)) +\
            '/../../data/protein-coding_gene.txt' 
try:
    fh = open(hgnc_file, 'rt')
    rd = csv.reader(fh, delimiter='\t')
    hgnc_dict = {}
    for row in rd:
        hgnc_id = row[0][5:]
        hgnc_name = row[1]
        hgnc_dict[hgnc_id] = hgnc_name
except IOError:
    hgnc_dict = {}

@lru_cache(maxsize=1000)
def get_hgnc_name(hgnc_id):
    try:
        hgnc_name = hgnc_dict[hgnc_id]
    except KeyError:
        xml_tree = get_hgnc_entry(hgnc_id)
        if xml_tree is None:
            return None
        hgnc_name_tag =\
            xml_tree.find("result/doc/str[@name='symbol']")
        if hgnc_name_tag is None:
            return None
        hgnc_name = hgnc_name_tag.text.strip()
    return hgnc_name

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
