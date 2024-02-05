from __future__ import print_function, unicode_literals
import logging
import os
import sys
__version__ = '1.22.0'

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

logger = logging.getLogger('indra')

from .config import get_config, has_config
