"""
This module provides an input interface and processor to the ISI reading
system.

The reader is set up to run within a Docker container.
For the ISI reader to run, set the Docker memory and swap space to the maximum.
For processing nxml files, install the nxml2txt utility
(https://github.com/spyysalo/nxml2txt) and set the configuration variable
NXML2TXT_PATH to its location. In addition, since the reader works with
Python 2 only, make sure PYTHON2_PATH is set in your config file or
environment and points to a Python 2 executable.
"""

from .api import process_text, process_nxml, process_preprocessed, \
    process_output_folder, process_json_file
