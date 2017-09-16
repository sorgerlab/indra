from indra.sources import trips
from indra.assemblers import PybelAssembler
import pybel
import requests
from indra.databases import hgnc_client
from indra.statements import *

def gnd(gene_name):
    return {'HGNC': hgnc_client.get_hgnc_id(gene_name)}

def ag(gene_name, **kwargs):
    return Agent(gene_name, db_refs=gnd(gene_name), **kwargs)

sos1 = ag('SOS1')
rasa1 = ag('RASA1')
kras = ag('KRAS')
braf = ag('BRAF', activity=ActivityCondition('kinase', True))
map2k1 = ag('MAP2K1')
map2k1_p = ag('MAP2K1', mods=[ModCondition('phosphorylation', 'S', '218')])
mapk1 = ag('MAPK1')
mapk1_p = ag('MAPK1', mods=[ModCondition('phosphorylation', 'T', '185')])
elk1 = ag('ELK1')
elk1_p = ag('ELK1', mods=[ModCondition('phosphorylation')])
elk1_tscript = ag('ELK1', activity=ActivityCondition('transcription', True))
fos = ag('FOS')

stmts = [
    Gef(sos1, kras),
    Gap(rasa1, kras),
    GtpActivation(kras, braf, 'kinase'),
    Phosphorylation(braf, map2k1, 'S', '218'),
    Phosphorylation(map2k1_p, mapk1, 'T', '185'),
    Phosphorylation(mapk1_p, elk1),
    ActiveForm(elk1_p, 'transcription', True),
    IncreaseAmount(elk1_tscript, fos),
]
#model_text = """
#    BRAF phosphorylates MAP2K1 at S218.
#    MAP2K1 phosphorylated at S218 phosphorylates MAPK1 at T185.
#    MAPK1 phosphorylated at T185 phosphorylates ELK1.
#    Phosphorylated ELK1 transcribes c-FOS.
#"""
#model_description = ' '.join([line.strip() for line in model_text.split('\n')])
#print("Processing text with TRIPS...")
#tp = trips.process_text(model_text)
model_description = 'Test of INDRA Statement assembly into PyBEL.'
print("Assembling to PyBEL...")

pba = PybelAssembler(stmts, name='INDRA_PyBEL_test',
                     description=model_description, version='0.0.7')
pba.make_model()
with open('pybel_model.json', 'wt') as f:
    pybel.to_json_file(pba.model, f)
url =  'https://pybel.scai.fraunhofer.de/api/receive'
headers = {'content-type': 'application/json'}
requests.post(url, json=pybel.to_json(pba.model), headers=headers)
