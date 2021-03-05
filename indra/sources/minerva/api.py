from .processor import SifProcessor



def process_file(filename, model_id):
    return process_files({model_id: filename})


def process_files(ids_to_filenames):
    model_id_to_sif_strs = {}
    for model_id, filename in ids_to_filenames.items():
        with open(filename, 'r') as f:
            sif_strs = f.readlines()
        model_id_to_sif_strs[model_id] = sif_strs
    return process_sif_strs(model_id_to_sif_strs)


def process_from_web(url):
    pass


def process_sif_strs(model_id_to_sif_strs):
    sp = SifProcessor(model_id_to_sif_strs)
    sp.extract_statements()
    return sp
