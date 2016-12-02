from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging

logger = logging.getLogger('sparser')


class SparserProcessor(object):
    def __init__(self, xml_etree):
        self.tree = xml_etree
        self.statements = []

