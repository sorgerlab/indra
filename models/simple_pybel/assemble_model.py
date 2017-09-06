from indra.sources import trips
from indra.assemblers import PybelAssembler
import pybel
import requests

model_text = """
    BRAF phosphorylates MAP2K1 at S218.
    MAP2K1 phosphorylated at S218 phosphorylates MAPK1 at T185.
    MAPK1 phosphorylated at T185 phosphorylates ELK1.
    Phosphorylated ELK1 transcribes c-FOS.
"""
model_description = ' '.join([line.strip() for line in model_text.split('\n')])
print("Processing text with TRIPS...")
tp = trips.process_text(model_text)
print("Assembling to PyBEL...")
pba = PybelAssembler(tp.statements, name='INDRA_PyBEL_test',
                     description=model_description, version='0.0.6')
pba.make_model()
with open('pybel_model.json', 'wt') as f:
    pybel.to_json_file(pba.model, f)
url =  'https://pybel.scai.fraunhofer.de/api/receive'
headers = {'content-type': 'application/json'}
requests.post(url, json=pybel.to_json(pba.model), headers=headers)

