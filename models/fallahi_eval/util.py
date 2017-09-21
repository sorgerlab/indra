from indra.util import _require_python3
import os
import json

# CREATE A JSON FILE WITH THIS INFORMATION, E.G., a file consisting of:
# {"basename": "fallahi_eval", "basedir": "output"}
with open('config.json', 'rt') as f:
    config = json.load(f)
# This is the base name used for all files created/saved
basen = config['basename']
# This is the base folder to read/write (potentially large) files from/to
# MODIFY ACCORDING TO YOUR OWN SETUP
based = config['basedir']

# This makes it easier to make standardized pickle file paths
prefixed_pkl = lambda suffix: os.path.join(based, basen + '_' + suffix + '.pkl')

