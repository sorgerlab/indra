from __future__ import print_function, unicode_literals
import logging
import os
import sys
__version__ = '1.15.0'

__all__ = ['assemblers', 'belief', 'databases', 'explanation', 'literature',
           'mechlinker', 'preassembler', 'sources', 'tools', 'util']

logging.basicConfig(format=('%(levelname)s: [%(asctime)s] %(name)s'
                            ' - %(message)s'),
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

from .config import get_config, has_config
