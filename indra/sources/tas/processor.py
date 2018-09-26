from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = []


class TasProcessor(object):
    """A processor for Target Affinity Spectrum data compiled by N. Moret.

    This data was compiled in the HMS LSP as an improvement on the "arbitrary"
    selection of targets present in the similar LINCS dataset.
    """
    def __init__(self, data):
        self._data = data
        self.statements = []
        for row in data:
            self._process_row(row)
        return

    def _process_row(self, row):
        pass
