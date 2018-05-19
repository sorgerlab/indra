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
hk1 = ag('HK1')
glu = Agent('D-GLUCOSE', db_refs={'CHEBI': 'CHEBI:17634'})
g6p = Agent('GLUCOSE-6-PHOSPHATE', db_refs={'CHEBI': 'CHEBI:4170'})
egfr = ag('EGFR')
grb2 = ag('GRB2')
tab1 = ag('TAB1')
p38_tab1 = Agent('MAPK14', bound_conditions=[BoundCondition(tab1)],
                 db_refs=gnd('MAPK14'))
egfr_dimer = Agent('EGFR', bound_conditions=[BoundCondition(egfr)],
                   db_refs=gnd('EGFR'))

sos1_bound = Agent('SOS1',
                   bound_conditions=[BoundCondition(egfr),
                                     BoundCondition(grb2)],
                   db_refs=gnd('SOS1'))
stmts = [
    Gef(sos1, kras),
    Gap(rasa1, kras),
    GtpActivation(kras, braf, 'kinase'),
    Phosphorylation(braf, map2k1, 'S', '218'),
    Phosphorylation(map2k1_p, mapk1, 'T', '185'),
    Phosphorylation(mapk1_p, elk1),
    ActiveForm(elk1_p, 'transcription', True),
    IncreaseAmount(elk1_tscript, fos),
    Conversion(hk1, [glu], [g6p]),
    Complex([egfr, grb2, sos1]),
    ActiveForm(sos1_bound, 'gef', True),
    Autophosphorylation(p38_tab1, 'Y', '100'),
    Transphosphorylation(egfr_dimer, 'Y', '1173'),
]

ev = Evidence('assertion', 'assertion', 'assertion', 'assertion')
for stmt in stmts:
    stmt.evidence = [ev]

model_description = 'Test of INDRA Statement assembly into PyBEL.'
print("Assembling to PyBEL...")

pba = PybelAssembler(stmts, name='INDRA_PyBEL_test',
                     description=model_description,
                     authors='John Bachman',
                     )
belgraph = pba.make_model()

# Write to BEL file
#pybel.to_bel_path(belgraph, 'simple_pybel.bel')

# Upload to PyBEL web
#with open('pybel_model.json', 'wt') as f:
#    pybel.to_json_file(pba.model, f)
#url =  'https://pybel.scai.fraunhofer.de/api/receive'
#headers = {'content-type': 'application/json'}
#requests.post(url, json=pybel.to_json(pba.model), headers=headers)

# Put in local database
pybel.to_database(pba.model)

