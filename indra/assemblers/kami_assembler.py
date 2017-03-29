from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

class KamiAssembler(object):
    def __init__(self, stmts=None):
        pass

    def add_statements(self, stmts):
        for stmt in stmts:
            self.statements.append(stmt)

