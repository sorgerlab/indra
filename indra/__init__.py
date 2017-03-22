from __future__ import print_function, unicode_literals
import logging
__version__ = '1.4.2'

__all__ = ['bel', 'biopax', 'trips', 'reach', 'index_cards', 'sparser',
           'databases', 'literature',
           'preassembler', 'assemblers', 'mechlinker', 'belief',
           'tools', 'util']
'''
#############
# For now these imports are disabled because
# (1) Every import would load everything in INDRA which is time consuming and
# (2) Optional dependencies in some modules will try to be loaded even if
# they are not intended to be used
##################
# Core
import statements
# Input processors
from indra import bel
from indra import biopax
from indra import trips
from indra import reach
from indra import index_cards
# Clients
from indra import databases
from indra import literature
# Assemblers
from indra import preassembler
from indra import assemblers
from indra import mechlinker
from indra import belief
# Tools and utils
from indra import tools
from indra import util
'''

logging.basicConfig(format='%(levelname)s: indra/%(name)s - %(message)s',
                    level=logging.INFO)

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
