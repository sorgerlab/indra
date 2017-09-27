import os
import shutil

__all__ = [
    'copy_default_config',
]

HERE = os.path.dirname(os.path.realpath(__file__))


def copy_default_config(destination):
    """Copies the default configuration to the given destination

    Parameters
    ----------
    destination : str
        The location to which a default RAS Machine config file is placed.
    """
    shutil.copy(os.path.join(HERE, 'default-config.yaml'), destination)
