import os
from pathlib import Path

default_bucket = 'world-modelers'
default_key_base = 'indra_models'
file_defaults = {'raw': 'raw_statements',
                 'sts': 'statements',
                 'cur': 'curations',
                 'meta': 'metadata'}


default_profile = 'wm'
HERE = Path(os.path.abspath(__file__)).parent
CACHE = HERE.joinpath('_local_cache')
CACHE.mkdir(exist_ok=True)


class InvalidCorpusError(Exception):
    pass


