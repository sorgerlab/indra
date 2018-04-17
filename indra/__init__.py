from __future__ import print_function, unicode_literals
import logging
import os
import shutil
import sys
if sys.version_info[0] == 3:
    from configparser import RawConfigParser
else:
    from ConfigParser import RawConfigParser
__version__ = '1.7.0'

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
config_path = os.path.join(config_dir, 'config.ini')
default_config_path = os.path.join(os.path.dirname(__file__),
                                   'resources/default_config.ini')
if not os.path.isfile(config_path):
    try:
        os.mkdir(config_dir)
    except Exception:
        logger.warning(config_dir + ' already exists')
    try:
        shutil.copyfile(default_config_path, config_path)
    except Exception:
        logger.warning('Could not copy default config file.')

# Load the configuration file into the config_file dictionary
# A ConfigParser-style configuration file can have multiple sections
# We ignore the section distinction  and load the key/value pairs from all
# sectionts into a single key/value list.

# Load key/value pairs from all sections into this dictionary
config_file = {}

try:
    parser = RawConfigParser()
    parser.optionxform = lambda x: x
    parser.read(config_path)
    sections = parser.sections()
    for section in sections:
        options = parser.options(section)
        for option in options:
            config_file[option] = str(parser.get(section, option))
except Exception:
    logger.warning('Could not load configuration file, will check ' +
                   'environment variables only for configuration.')
    config_file = {}

# Expand ~ to the home directory
for key in config_file:
    config_file[key] = os.path.expanduser(config_file[key])

# Some logic in INDRA checks to see if the config file is None if
# omitted. Set config values to None if the empty string.
for key in config_file:
    if config_file[key] == "":
        config_file[key] = None


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


def has_config(key):
    """Returns whether the configuration value for the given kehy is present.

    Parmeters
    ---------
    key: str
        The key for the configuration value to fetch

    Returns
    -------
    value: bool
        Whether the configuration value is present
    """
    return get_config(key) is not None
