import pickle
import tempfile
import os
import os.path
import subprocess

nxml_dir = '/Users/daniel/Downloads/nxml2txt'
xml_path = os.path.join(nxml_dir, 'nxml2txt')

xmls = pickle.load(open('xml_samples.p', 'rb'))

for xml in xmls:
    fname = tempfile.mkstemp('_indra_xml_to_text')

    with open(fname, 'w') as f:
        f.write(xml)

    os.remove(fname)


