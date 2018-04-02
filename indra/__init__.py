from __future__ import print_function, unicode_literals
import logging
import os
import json
import shutil
__version__ = '1.6.0'

__all__ = ['assemblers', 'belief', 'databases', 'explanation', 'literature',
           'mechlinker', 'preassembler', 'sources', 'tools', 'util']

logging.basicConfig(format='%(levelname)s: [%(asctime)s] indra/%(name)s - %(message)s',
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

# Suppress INFO-level logging from some dependencies
logging.getLogger('requests').setLevel(logging.ERROR)
logging.getLogger('urllib3').setLevel(logging.ERROR)
logging.getLogger('rdflib').setLevel(logging.ERROR)
logging.getLogger('boto3').setLevel(logging.CRITICAL)
logging.getLogger('botocore').setLevel(logging.CRITICAL)

# This is specifically to suppress lib2to3 logging from networkx
import lib2to3.pgen2.driver
class Lib2to3LoggingModuleShim(object):
    def getLogger(self):
        return logging.getLogger('lib2to3')
lib2to3.pgen2.driver.logging = Lib2to3LoggingModuleShim()
logging.getLogger('lib2to3').setLevel(logging.ERROR)

logger = logging.getLogger('indra')

# If the configuration file does not exist, try to create it from the default
config_dir = os.path.expanduser('~/.config/indra')
config_path = os.path.join(config_dir, 'config.json')
if not os.path.isfile(config_path):
    os.mkdir(config_dir)
    default_config = os.path.join(os.path.dirname(__file__),
                                  'resources/default_config.json')
    shutil.copyfile(default_config, config_path)

# Load the configuration file
with open(config_path, 'r') as f:
    config_contents = f.read()
    config_file = json.loads(config_contents)

    # Expand ~ to the home directory
    for key in config_file:
        config_file[key] = os.path.expanduser(config_file[key])

def get_config(key):
    """Returns the configuration value, first checking the environemnt
    variables and then, if it's not present there, checking the configuration
    file.

    Parameters
    ----------
    key: str
        The key for the configuration value to fetch

    Returns
    -------
    value: str
        The configuration value
    """
    if key in os.environ:
        return os.environ[key]
    elif key in config_file:
        return config_file[key]
    else:
        logger.warning('Could not find ' + str(key) +
                           ' in environment variables or configuration file')
        return None
