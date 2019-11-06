from indra.tools.reading.readers.core import EmptyReader

from indra.sources import trips


class TripsReader(EmptyReader):
    """A stand-in for TRIPS reading.

    Currently, we do not run TRIPS (more specifically DRUM) regularly at large
    scales, however on occasion we have outputs from TRIPS that were generated
    a while ago.
    """
    name = 'TRIPS'

    def __init__(self, *args, **kwargs):
        self.version = self.get_version()
        return

    def _read(self, *args, **kwargs):
        return []

    @classmethod
    def get_version(cls):
        return 'STATIC'

    @staticmethod
    def get_processor(content):
        return trips.process_xml(content)



