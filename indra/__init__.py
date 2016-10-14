from __future__ import print_function, unicode_literals
import logging
__version__ = '1.3.0'

logging.basicConfig(format='%(levelname)s: indra/%(name)s - %(message)s',
                    level=logging.INFO)

logging.getLogger('requests').setLevel(logging.ERROR)
logging.getLogger('urllib3').setLevel(logging.ERROR)
logging.getLogger('rdflib').setLevel(logging.ERROR)
logging.getLogger('boto3').setLevel(logging.CRITICAL)
logging.getLogger('botocore').setLevel(logging.CRITICAL)
