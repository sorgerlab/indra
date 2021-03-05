import requests
from .processor import SifProcessor


base_url = ('https://git-r3lab.uni.lu/covid/models/-/raw/master/'
            'Executable%20Modules/SBML_qual_build/sif/')
model_ids_to_filenames = {
    786: 'Virus_replication_cycle_stable_raw.sif',
    787: 'PAMP_signalling_stable_raw.sif',
    788: 'Interferon1_stable_raw.sif',
    789: 'Orf3a_stable_raw.sif',
    790: 'TGFB_pathway_stable_raw.sif',
    791: 'IFN-lambda_stable_raw.sif',
    792: 'Kynurenine_pathway_stable_raw.sif',
    793: 'HMOX1_Pathway_stable_raw.sif',
    794: 'Orf10_Cul2_pathway_stable_raw.sif',
    795: 'E_protein_stable_raw.sif',
    796: 'RTC-and-transcription_stable_raw.sif',
    797: 'JNK_pathway_stable_raw.sif',
    798: 'ER_Stress_stable_raw.sif',
    799: 'Apoptosis_stable_raw.sif',
    801: 'Coagulation-pathway_stable_raw.sif',
    802: 'Nsp4_Nsp6_stable_raw.sif',
    803: 'Pyrimidine_deprivation_stable_raw.sif',
    804: 'ETC_stable_raw.sif',
    805: 'Renin_angiotensin_stable_raw.sif',
    806: 'Nsp9_protein_stable_raw.sif',
    807: 'NLRP3_Activation_stable_raw.sif'
    }


def process_file(filename, model_id):
    return process_files({model_id: filename})


def process_files(ids_to_filenames):
    model_id_to_sif_strs = {}
    for model_id, filename in ids_to_filenames.items():
        with open(filename, 'r') as f:
            sif_strs = f.readlines()
        model_id_to_sif_strs[model_id] = sif_strs
    return process_sif_strs(model_id_to_sif_strs)


def process_from_web(model_ids='all'):
    if model_ids == 'all':
        model_ids = list(model_ids_to_filenames.keys())
    model_id_to_sif_strs = {}
    for model_id in model_ids:
        fname = model_ids_to_filenames[model_id]
        url = base_url + fname
        res = requests.get(url)
        if res.status_code == 200:
            sif_strs = res.text.split('\n')
            model_id_to_sif_strs[model_id] = sif_strs
    return process_sif_strs(model_id_to_sif_strs)


def process_sif_strs(model_id_to_sif_strs):
    sp = SifProcessor(model_id_to_sif_strs)
    sp.extract_statements()
    return sp
