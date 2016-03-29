import os
import json
import tempfile
import requests
from indra.java_vm import autoclass, JavaException
import indra.databases.pmc_client as pmc_client
import indra.databases.pubmed_client as pubmed_client
from processor import ReachProcessor
from reach_reader import ReachReader

reach_text_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/text'
reach_nxml_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/nxml'

# For offline reading
reach_reader = ReachReader()

def process_pmc(pmc_id, offline=False):
    xml_str = pmc_client.get_xml(pmc_id)
    if xml_str is None:
        return None
    fname = pmc_id + '.nxml'
    with open(fname, 'wt') as fh:
        fh.write(xml_str.encode('utf-8'))
    rp = process_nxml_str(xml_str, citation=pmc_id, offline=offline)
    return rp

def process_pubmed_abstract(pubmed_id, offline=False):
    abs_txt = pubmed_client.get_abstract(pubmed_id)
    if abs_txt is None:
        return None
    rp = process_text(abs_txt, citation=pubmed_id, offline=offline)
    return rp

def process_text(text, citation=None, offline=False):
    if offline:
        api_ruler = reach_reader.get_api_ruler()
        if api_ruler is None:
            print 'Cannot read offline because the REACH ApiRuler could not' +\
                    ' be instantiated.'
            return None
        try:
            result_map = api_ruler.annotateText(text, 'fries')
        except JavaException:
            print 'Could not process text.'
            return None
        json_str = result_map.get('resultJson')
    else:
        data = {'text': text.encode('utf-8')}
        try:
            res = requests.post(reach_text_url, data)
        except requests.exceptions.RequestException as e:
            print 'Could not connect to REACH service:'
            print e
            return None
        # TODO: we could use res.json() here to get a dict 
        # directly
        json_str = res.text
    with open('reach_output.json', 'wt') as fh:
        out_str = json_str
        fh.write(out_str)
    return process_json_str(json_str, citation)

def process_nxml_str(nxml_str, citation=None, offline=False):
    if offline:
        api_ruler = reach_reader.get_api_ruler()
        if api_ruler is None:
            print 'Cannot read offline because the REACH ApiRuler could not' +\
                    ' be instantiated.'
            return None
        try:
            #TODO: Test if UTF-8 files are parsed correctly here
            result_map = api_ruler.annotateNxml(nxml_str, 'fries')
        except JavaException:
            print 'Could not process NXML.'
            return None
        json_str = result_map.get('resultJson')
    else:
        data = {'nxml': nxml_str}
        try:
            res = requests.post(reach_nxml_url, data)
        except requests.exceptions.RequestException as e:
            print 'Could not connect to REACH service:'
            print e
            return None
        if res.status_code != 200:
            print 'Could not process NXML via REACH service.'
            return None
        json_str = res.text
    with open('reach_output.json', 'wt') as fh:
        fh.write(json_str)
    return process_json_str(json_str, citation)

def process_nxml_file(file_name, citation=None, offline=False):
    nxml_str = open(file_name, 'rt').read()
    nxml_str = nxml_str.decode('utf-8')
    return process_nxml_str(nxml_str, citation, offline)

def process_json_file(file_name, citation=None):
    try:
        with open(file_name, 'rt') as fh:
            json_str = fh.read()
            return process_json_str(json_str, citation)
    except IOError:
        print 'Could not read file %s.' % file_name


def process_json_str(json_str, citation=None):
    json_str = json_str.replace('frame-id','frame_id')
    json_str = json_str.replace('argument-label','argument_label')
    json_str = json_str.replace('object-meta','object_meta')
    json_str = json_str.replace('doc-id','doc_id')
    json_str = json_str.replace('is-hypothesis','is_hypothesis')
    json_str = json_str.replace('is-negated','is_negated')
    json_dict = json.loads(json_str)
    rp = ReachProcessor(json_dict, citation)
    rp.get_phosphorylation()
    rp.get_complexes()
    rp.get_activation()
    return rp

if __name__ == '__main__':
    rp = process_json_file('PMC0000001.uaz.events.json')
