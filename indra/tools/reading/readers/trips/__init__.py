import os
import logging
import subprocess as sp
from datetime import datetime, timedelta, timezone

from indra.tools.reading.readers.core import Reader

from indra.sources.trips import client, process_xml


logger = logging.getLogger(__name__)

startup_path = '/sw/drum/bin/startup.sh'
service_host = 'drum'
DRUM_DOCKER = '292075781285.dkr.ecr.us-east-1.amazonaws.com/drum'


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
            logger.info("Starting up a TRIPS service from within the docker.")
            p = sp.Popen([startup_path], stdout=sp.PIPE,
                         stderr=sp.STDOUT)
            service_endpoint = 'http://localhost:80/cgi/'
        else:
            logger.info("Starting up a TRIPS service using drum docker.")
            p = sp.Popen(['docker', 'run', '-it', '-p', '8080:80',
                          '--entrypoint', startup_path, DRUM_DOCKER],
                         stdout=sp.PIPE, stderr=sp.STDOUT)
            service_endpoint = 'http://localhost:8080/cgi/'

        # Wait for the service to be ready
        for line in iter(p.stdout.readline, b''):
            log_line = line.strip().decode('utf8')
            logger.info('TRIPS: ' + log_line)
            if log_line == 'ready.':
                break

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
        git_date_cmd = ['git', 'log', '-1', '--format=%cd']
        if os.environ.get("IN_TRIPS_DOCKER", "false") == "true":
            curdir = os.getcwd()
            try:
                # Get the date of the last commit from git in drum.
                os.chdir('/sw/drum')
                res = sp.run(git_date_cmd, stdout=sp.PIPE)
            finally:
                os.chdir(curdir)
        else:
            # We need to add a function within the docker to
            assert False, "This is currently not available."

        # Format that string into a datetime and standardize to utc.
        d = datetime.strptime(res.stdout.decode('utf-8').stirp(),
                              '%a %b %d %H:%M:%S %Y %z')
        d.astimezone(tz=timezone(timedelta(0)))

        # Create a formatted string as the version.
        version = d.strftime('%Y%b%d')

        return version

    @staticmethod
    def get_processor(reading_content):
        return process_xml(reading_content)
