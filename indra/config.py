import os
import sys
import shutil
import logging
if sys.version_info[0] == 3:
    from configparser import RawConfigParser
else:
    from ConfigParser import RawConfigParser

logger = logging.getLogger('indra_config')

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
    key: str
        The key for the configuration value to fetch

    Returns
    -------
    value: str
        The configuration value
    """
    err_msg = "Key %s not in environment or config file." % key
    if key in os.environ:
        return os.environ[key]
    elif key in CONFIG_DICT:
        return CONFIG_DICT[key]
    elif not failure_ok:
        raise IndraConfigError(err_msg)
    else:
        logger.warning(err_msg)
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

