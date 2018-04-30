import sys
import logging
from os import path, mkdir, environ
from shutil import copyfile
if sys.version_info[0] == 3:
    from configparser import ConfigParser
else:
    from ConfigParser import ConfigParser

FILE_PATH = path.dirname(path.abspath(__file__))
DEFAULT_DB_CONFIG_PATH = path.join(FILE_PATH, 'default_db_config.ini')

DB_CONFIG_DIR = path.expanduser('~/.config/indra')
DB_CONFIG_PATH = path.join(DB_CONFIG_DIR, 'db_config.ini')

DB_STR_FMT = "{prefix}://{username}:{password}{host}{port}/{name}"
ENV_PREFIX = 'INDRADB'


if not path.exists(DB_CONFIG_DIR):
    mkdir(DB_CONFIG_DIR)


if not path.exists(DB_CONFIG_PATH):
    copyfile(DEFAULT_DB_CONFIG_PATH, DB_CONFIG_PATH)

DATABASES = None


def get_databases(force_update=False):
    global DATABASES
    if DATABASES is None or force_update:
        DATABASES = {}
        parser = ConfigParser()
        parser.read(DB_CONFIG_PATH)
        for db_name in parser.sections():
            def_dict = {k: parser.get(db_name, k) for k in parser.options(db_name)}
            if def_dict['host']:
                def_dict['host'] = '@' + def_dict['host']
            def_dict['prefix'] = def_dict['dialect']
            if def_dict['driver']:
                def_dict['prefix'] += def_dict['driver']
            if def_dict['port']:
                def_dict['port'] += ':' + def_dict['port']
            DATABASES[db_name] = DB_STR_FMT.format(**def_dict)
        DATABASES.update({k[len(ENV_PREFIX):].lower(): v
                          for k, v in environ.items()
                          if k.startswith(ENV_PREFIX)})
    return DATABASES
