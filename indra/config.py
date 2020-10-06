import os
import sys
import shutil
import logging
if sys.version_info[0] == 3:
    from configparser import RawConfigParser
else:
    from ConfigParser import RawConfigParser

logger = logging.getLogger(__name__)

# If the configuration file does not exist, try to create it from the default
home_dir = os.path.expanduser('~')
config_dir = os.path.join(home_dir, '.config', 'indra')
config_path = os.path.join(config_dir, 'config.ini')
default_config_path = os.path.join(os.path.dirname(__file__),
                                   'resources/default_config.ini')
if not os.path.isfile(config_path):
    try:
        os.makedirs(config_dir)
    except Exception:
        logger.warning(config_dir + ' already exists')
    try:
        shutil.copyfile(default_config_path, config_path)
    except Exception:
        logger.warning('Could not copy default config file.')


def populate_config_dict(config_path):
    """Load the configuration file into the config_file dictionary

    A ConfigParser-style configuration file can have multiple sections, but
    we ignore the section distinction  and load the key/value pairs from all
    sections into a single key/value list.
    """
    try:
        config_dict = {}
        parser = RawConfigParser()
        parser.optionxform = lambda x: x
        parser.read(config_path)
        sections = parser.sections()
        for section in sections:
            options = parser.options(section)
            for option in options:
                config_dict[option] = str(parser.get(section, option))
    except Exception as e:
        logger.warning("Could not load configuration file due to exception. "
                       "Only environment variable equivalents will be used.")
        return None

    for key in config_dict.keys():
        if config_dict[key] == '':
            config_dict[key] = None
        elif isinstance(config_dict[key], str):
            config_dict[key] = os.path.expanduser(config_dict[key])
    return config_dict


def _check_config_dict():
    # Check the keys against the default.
    default_config_dict = populate_config_dict(default_config_path)
    for key in default_config_dict.keys():
        if key not in CONFIG_DICT:
            logger.debug("Key %s found in default config but not in %s."
                         % (key, config_path))
    for key in CONFIG_DICT.keys():
        if key not in default_config_dict:
            logger.debug("Key %s found in %s but not in default config."
                         % (key, config_path))


CONFIG_DICT = populate_config_dict(config_path)
if CONFIG_DICT is None:
    CONFIG_DICT = {}
else:
    _check_config_dict()


class IndraConfigError(Exception):
    pass


def get_config(key, failure_ok=True):
    """Get value by key from config file or environment.

    Returns the configuration value, first checking the environment
    variables and then, if it's not present there, checking the configuration
    file.

    Parameters
    ----------
    key : str
        The key for the configuration value to fetch
    failure_ok : Optional[bool]
        If False and the configuration is missing, an IndraConfigError is
        raised. If True, None is returned and no error is raised in case
        of a missing configuration. Default: True

    Returns
    -------
    value : str or None
        The configuration value or None if the configuration value doesn't
        exist and failure_ok is set to True.
    """
    err_msg = "Key %s not in environment or config file." % key
    if key in os.environ:
        return os.environ[key]
    elif key in CONFIG_DICT:
        val = CONFIG_DICT[key]
        # We interpret an empty value in the config file as a failure
        if val is None and not failure_ok:
            msg = 'Key %s is set to an empty value in config file.' % key
            raise IndraConfigError(msg)
        else:
            return val
    elif not failure_ok:
        raise IndraConfigError(err_msg)
    else:
        return None


def has_config(key):
    """Returns whether the configuration value for the given kehy is present.

    Parameters
    ----------
    key : str
        The key for the configuration value to fetch

    Returns
    -------
    value : bool
        Whether the configuration value is present
    """
    return get_config(key, failure_ok=True) is not None


KNOWLEDGE_SOURCE_INFO = {
    'eidos': {
        'name': 'Eidos',
        'link': 'https://github.com/clulab/eidos',
        'type': 'reader',
        'domain': 'general'
    },
    'cwms': {
        'name': 'TRIPS/CWMS',
        'link': 'http://trips.ihmc.us/parser/cgi/cwmsreader',
        'type': 'reader',
        'domain': 'general'
    },
    'hume': {
        'name': 'Hume',
        'link': 'https://github.com/BBN-E/Hume',
        'type': 'reader',
        'domain': 'general'
    },
    'sofia': {
        'name': 'Sofia',
        'link': 'https://sofia.worldmodelers.com/ui/',
        'type': 'reader',
        'domain': 'general'
    },
    'drum': {
        'name': 'TRIPS/DRUM',
        'link': 'http://trips.ihmc.us/parser/cgi/drum',
        'type': 'reader',
        'domain': 'biology'
    },
    'reach': {
        'name': 'REACH',
        'link': 'https://github.com/clulab/reach',
        'type': 'reader',
        'domain': 'biology'
    },
    'sparser': {
        'name': 'Sparser',
        'link': 'https://github.com/ddmcdonald/sparser',
        'type': 'reader',
        'domain': 'biology'
    },
    'tees': {
        'name': 'TEES',
        'link': 'https://github.com/jbjorne/TEES',
        'type': 'reader',
        'domain': 'biology'
    },
    'medscan': {
        'name': 'MedScan',
        'link': 'https://doi.org/10.1093/bioinformatics/btg207',
        'type': 'reader',
        'domain': 'biology'
    },
    'rlimsp': {
        'name': 'RLIMS-P',
        'link': 'https://research.bioinformatics.udel.edu/rlimsp',
        'type': 'reader',
        'domain': 'biology'
    },
    'isi': {
        'name': 'ISI/AMR',
        'link': 'https://github.com/sgarg87/big_mech_isi_gg',
        'type': 'reader',
        'domain': 'biology'
    },
    'geneways': {
        'name': 'Geneways',
        'link': 'https://www.ncbi.nlm.nih.gov/pubmed/15016385',
        'type': 'reader',
        'domain': 'biology'
    },
    'biopax': {
        'name': 'BioPax',
        'link': 'http://www.biopax.org/',
        'type': 'database',
        'domain': 'biology'
    },
    'pc': {
        'name': 'PathwayCommons',
        'link': 'http://pathwaycommons.org/',
        'type': 'database',
        'domain': 'biology'
    },
    'bel': {
        'name': 'BEL',
        'link': 'https://github.com/pybel/pybel',
        'type': 'database',
        'domain': 'biology'
    },
    'signor': {
        'name': 'Signor',
        'link': 'https://signor.uniroma2.it/',
        'type': 'database',
        'domain': 'biology'
    },
    'biogrid': {
        'name': 'BioGRID',
        'link': 'https://thebiogrid.org/',
        'type': 'database',
        'domain': 'biology'
    },
    'tas': {
        'name': 'Target Affinity Spectrum',
        'link': 'https://doi.org/10.1101/358978',
        'type': 'database',
        'domain': 'biology'
    },
    'lincs_drug': {
        'name': 'LINCS small molecules',
        'link': 'http://lincs.hms.harvard.edu/db/sm/',
        'type': 'database',
        'domain': 'biology'
    },
    'hprd': {
        'name': 'HPRD',
        'link': 'http://www.hprd.org',
        'type': 'database',
        'domain': 'biology'
    },
    'trrust': {
        'name': 'TRRUST',
        'link': 'https://www.grnpedia.org/trrust/',
        'type': 'database',
        'domain': 'biology'
    },
    'phosphoelm': {
        'name': 'Phospho.ELM',
        'link': 'http://phospho.elm.eu.org/',
        'type': 'database',
        'domain': 'biology'
    },
    'virhostnet': {
        'name': 'VirHostNet',
        'link': 'http://virhostnet.prabi.fr/',
        'type': 'database',
        'domain': 'biology'
    },
    'ctd': {
        'name': 'CTD',
        'link': 'http://ctdbase.org',
        'type': 'database',
        'domain': 'biology'
    },
    'drugbank': {
        'name': 'DrugBank',
        'link': 'https://www.drugbank.ca/',
        'type': 'database',
        'domain': 'biology'
    },
    'omnipath': {
        'name': 'OmniPath',
        'link': 'https://omnipathdb.org/',
        'type': 'database',
        'domain': 'biology'
    }
}
