import os
import subprocess as sp

from indra.tools.reading.readers.core import Reader

from indra.sources.trips import client, process_xml
from indra_db import formats


class TripsReader(Reader):
    """A stand-in for TRIPS reading.

    Currently, we do not run TRIPS (more specifically DRUM) regularly at large
    scales, however on occasion we have outputs from TRIPS that were generated
    a while ago.
    """
    name = 'TRIPS'
    result_format = formats.XML

    def __init__(self, *args, **kwargs):
        self.version = self.get_version()
        return

    def _read(self, content_iter, verbose=False, log=False, n_per_proc=None):
        # Start trips running
        if os.environ.get("IN_TRIPS_DOCKER", 'false') != 'true':
            return []

        p = sp.Popen('/sw/drum/bin/startup.sh', stdout=sp.PIPE,
                     stderr=sp.STDOUT)
        service_endpoint = 'http://localhost:80/cgi/'
        service_host = 'drum'

        # Process all the content.
        for content in content_iter:
            html = client.send_query(content.get_text(),
                                     service_endpoint=service_endpoint,
                                     service_host=service_host)
            xml = client.get_xml(html)
            self.add_result(content.get_id(), xml)

        return self.results

    @classmethod
    def get_version(cls):
        return 'STATIC'

    @staticmethod
    def get_processor(reading_content):
        return process_xml(reading_content)



