import os
import json
import tempfile
import urllib, urllib2
import requests
from indra.java_vm import autoclass, JavaException
import indra.databases.pmc_client as pmc_client
import indra.databases.pubmed_client as pubmed_client
from processor import ReachProcessor

reach_text_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/text'
reach_nxml_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/nxml'

def process_pmc(pmc_id):
    if pmc_id.upper().startswith('PMC'):
        pmc_id = pmc_id[3:]
    xml_str = pmc_client.get_xml(pmc_id)
    if xml_str is None:
        return None
    with tempfile.NamedTemporaryFile() as fh:
        fh.write(xml_str)
        fh.flush()
        rp = process_nxml(fh.name)
    return rp

def process_pubmed_abstract(pubmed_id, offline=False):
    if pubmed_id.upper().startswith('PMID'):
        pubmed_id = pubmed_id[4:]
    abs_txt = pubmed_client.get_abstract(pubmed_id)
    if abs_txt is None:
        return None
    rp = process_text(abs_txt, offline)
    return rp

def process_text(txt, offline=False):
    if offline:
        nxml_txt = '<article><body><sec><p>%s</p></sec></body></article>' % txt
        tmp_file = tempfile.NamedTemporaryFile()
        tmp_file.file.write(nxml_txt)
        tmp_file.file.flush()
        return process_nxml(tmp_file.name)
    else:
        req = urllib2.Request(reach_text_url, data=urllib.urlencode({'text': txt}))
        res = urllib2.urlopen(req)
        json_str = res.read()
        #json_dict = json.loads(json_str)
        #events_dict = json_dict['events']
        #events_json_str = json.dumps(events_dict, indent=1)
        with open('reach_output.json', 'wt') as fh:
            fh.write(json_str)
        return process_json_str(json_str)

def process_nxml(file_name, offline=False):
    if offline:
        try:
            api_ruler = autoclass('edu.arizona.sista.reach.apis.ApiRuler')
            result_map = api_ruler.annotateNxml(file_name, 'fries')
        except JavaException:
            print 'Could not process file %s.' % file_name
            return None
        json_str = result_map.get('resultJson')
    else:
        txt = open(file_name, 'rt').read()
        req = urllib2.Request(reach_nxml_url, data=urllib.urlencode({'nxml': txt}))
        res = urllib2.urlopen(req)
        json_str = res.read()
    with open('reach_output.json', 'wt') as fh:
        fh.write(json_str)
    return process_json_str(json_str)


def process_json_file(file_name):
    try:
        with open(file_name, 'rt') as fh:
            json_str = fh.read()
            return process_json_str(json_str)
    except IOError:
        print 'Could not read file %s.' % file_name


def process_json_str(json_str):
    json_str = json_str.replace('frame-id','frame_id')
    json_str = json_str.replace('argument-label','argument_label')
    json_str = json_str.replace('object-meta','object_meta')
    json_str = json_str.replace('doc-id','doc_id')
    json_dict = json.loads(json_str)
    rp = ReachProcessor(json_dict)
    rp.get_phosphorylation()
    rp.get_complexes()
    rp.get_activation()
    return rp

if __name__ == '__main__':
    rp = process_json_file('PMC0000001.uaz.events.json')
