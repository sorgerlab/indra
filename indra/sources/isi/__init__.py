"""
This module provides an input interface and processor to the ISI reading
system.

The reader is set up to run within a Docker container.
For the ISI reader to run, set the Docker memory and swap space to the maximum.
"""

from .api import process_text, process_nxml, process_preprocessed, \
    process_output_folder, process_json_file
