"""
The BioOntology is built and cached automatically during runtime. If
a cached version already exists, it is loaded from the cache.

To control the build and clean up caches if necessary, one can call

.. code-block:: bash

    python -m indra.ontology.bio <operation>

to build or clean up the INDRA bio ontology. The script takes
a single operation argument which can be as follows:

* `build`: build the ontology and cache it
* `clean`: delete the current version of the ontology from the cache
* `clean-old`: delete all versions of the ontology except the current one
* `clean-all`: delete all versions of the bio ontology from the cache
"""
import os
import sys
import glob
import shutil
import logging
from .ontology import BioOntology, CACHE_DIR

logger = logging.getLogger('indra.ontology.bio')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        logger.info('Operation missing. Supported operations: '
                    'build, clean, clean-old, clean-all.')
        sys.exit(1)
    operation = sys.argv[1]
    if operation == 'build':
        BioOntology().initialize()
    elif operation.startswith('clean'):
        parent_dir = os.path.normpath(os.path.join(CACHE_DIR, os.pardir))
        version_paths = glob.glob(os.path.join(parent_dir, '*', ''))
        if operation == 'clean-all':
            to_remove = [parent_dir]
        else:
            to_remove = []
            for version_path in version_paths:
                version = os.path.basename(os.path.normpath(version_path))
                if operation == 'clean-old' and version != BioOntology.version:
                    to_remove.append(version_path)
                elif operation == 'clean' and version == BioOntology.version:
                    to_remove.append(version_path)
        for rem in to_remove:
            logger.info('Removing %s' % rem)
            shutil.rmtree(rem)
