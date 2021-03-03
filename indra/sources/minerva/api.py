from .processor import SifProcessor


def process_file(filename):
    with open(filename, 'r') as f:
        sif_strs = f.readlines()
    return process_sif_strs(sif_strs)


def process_files(filenames):
    sif_strs = []
    for filename in filenames:
        with open(filename, 'r') as f:
            strs = f.readlines()
        sif_strs += strs
    return process_sif_strs(sif_strs)


def process_from_web(url):
    pass


def process_sif_strs(sif_strs):
    sp = SifProcessor(sif_strs)
    sp.extract_statements()
    return sp
