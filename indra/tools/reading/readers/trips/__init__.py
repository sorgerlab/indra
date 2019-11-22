import os
import subprocess as sp

from indra.tools.reading.readers.core import Reader

from indra.sources.trips import client, process_xml

startup_path = '/sw/drum/bin/startup.sh'
service_host = 'drum'


class TripsReader(Reader):
    """A stand-in for TRIPS reading.

    Currently, we do not run TRIPS (more specifically DRUM) regularly at large
    scales, however on occasion we have outputs from TRIPS that were generated
    a while ago.
    """
    name = 'TRIPS'
    result_format = 'xml'

    def __init__(self, *args, **kwargs):
        self.version = self.get_version()
        super(TripsReader, self).__init__(*args, **kwargs)
        return

    def _read(self, content_iter, verbose=False, log=False, n_per_proc=None):
        # Start trips running
        if os.environ.get("IN_TRIPS_DOCKER", 'false') == 'true':
            p = sp.Popen([startup_path], stdout=sp.PIPE,
                         stderr=sp.STDOUT)
            service_endpoint = 'http://localhost:80/cgi/'
        else:
            p = sp.Popen(['docker', 'run', '-id', '-p', '8080:80',
                          '--entrypoint', startup_path, 'drum:latest'])
            service_endpoint = 'http://localhost:8080/cgi/'

        # Process all the content.
        for content in content_iter:
            html = client.send_query(content.get_text(),
                                     service_endpoint=service_endpoint,
                                     service_host=service_host)
            xml = client.get_xml(html)
            self.add_result(content.get_id(), xml)

        # Stop TRIPS
        p.kill()

        return self.results

    @classmethod
    def get_version(cls):
        return 'STATIC'

    @staticmethod
    def get_processor(reading_content):
        return process_xml(reading_content)



