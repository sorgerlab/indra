import os
import socket
import random
import logging
import threading
import subprocess as sp
from unidecode import unidecode
from contextlib import closing
from datetime import datetime, timedelta, timezone

from indra.resources import greek_alphabet
from indra.tools.reading.readers.core import Reader

from indra.sources.trips import client, process_xml


logger = logging.getLogger(__name__)

startup_path = '/sw/drum/bin/startup.sh'
service_endpoint = 'drum'
DRUM_DOCKER = '292075781285.dkr.ecr.us-east-1.amazonaws.com/drum'


def find_free_ports():
    """Find ports that are unused.

    The order is randomized to minimize the chances of race-condition overlaps.
    """
    ports = list(range(1, 65536))
    random.shuffle(ports)
    for port in ports:
        with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sok:
            res = sok.connect_ex(('localhost', port))
            if res != 0:
                yield port


class TripsStartupError(Exception):
    pass


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
        self.running = False
        self.stopping = False
        return

    def _monitor_trips_service(self, proc):
        self.running = True
        for _ in _tail_trips(proc):
            if self.stopping:
                break
        proc.communicate()
        if proc.returncode:
            logger.error("TRIPS ended badly.")
        self.running = False

    def _read(self, content_iter, verbose=False, log=False, n_per_proc=None):
        # Start trips running
        for port in find_free_ports():
            port_failure = False
            if os.environ.get("IN_TRIPS_DOCKER", 'false') == 'true':
                logger.info("Attempting to starting up a TRIPS service from "
                            "within the docker on port %d." % port)
                p = sp.Popen([startup_path], stdout=sp.PIPE,
                             stderr=sp.STDOUT)
                service_host = 'http://localhost:%d/cgi/' % port
            else:
                logger.info("Starting up a TRIPS service using drum docker.")
                p = sp.Popen(['docker', 'run', '-it', '-p', '8080:%d' % port,
                              '--entrypoint', startup_path, DRUM_DOCKER],
                             stdout=sp.PIPE, stderr=sp.STDOUT)
                service_host = 'http://localhost:8080/cgi/'

            # Wait for the service to be ready
            for log_line in _tail_trips(p):
                if 'can\'t bind to port' in 'log_line':
                    port_failure = True
                if log_line == 'Ready':
                    # TRIPS is ready to read and we can continue on.
                    break
            else:
                logger.error("TRIPS failed to start on port %d." % port)
                if port_failure:
                    # If the failure due to wrong port, try another port.
                    continue
                else:
                    # Otherwise give up.
                    raise TripsStartupError("Trips failed to start up.")

            # The above for-loop existed without going to else, indicating
            # successful startup.
            break
        else:
            # We exhausted all possible ports.
            raise TripsStartupError("Could not start TRIPS, all ports appear "
                                    "to be busy.")
        logger.info("Service has started up.")

        # Set up the trips monitor
        th = threading.Thread(target=self._monitor_trips_service, args=[p])
        th.start()

        # Process all the content.
        for content in content_iter:
            if not self.running:
                logger.error("Breaking loop: trips is down.")
                break

            # Clean up the text string a bit.
            raw_text = content.get_text()
            for greek_letter, spelled_letter in greek_alphabet.items():
                raw_text = raw_text.replace(greek_letter, spelled_letter)
            text = unidecode(raw_text)

            # Process the text
            html = client.send_query(text, service_host=service_host,
                                     service_endpoint=service_endpoint)
            xml = client.get_xml(html)
            self.add_result(content.get_id(), xml)

        # Stop TRIPS if it hasn't stopped already.
        p.kill()  # Send signal to the process to stop
        self.stopping = True  # Sends signal to the loop to stop
        th.join()

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
            res = sp.run(['docker', 'run', DRUM_DOCKER, '-c',
                          ' '.join(['cd', '/sw/drum;'] + git_date_cmd)],
                         stdout=sp.PIPE)

        # Format that string into a datetime and standardize to utc.
        d = datetime.strptime(res.stdout.decode('utf-8').strip(),
                              '%a %b %d %H:%M:%S %Y %z')
        d.astimezone(tz=timezone(timedelta(0)))

        # Create a formatted string as the version.
        version = d.strftime('%Y%b%d')

        return version

    @staticmethod
    def get_processor(reading_content):
        return process_xml(reading_content)


def _tail_trips(proc):
    for line in iter(proc.stdout.readline, b''):
        log_line = line.strip().decode('utf-8')
        logger.info('TRIPS: ' + log_line)
        yield log_line
