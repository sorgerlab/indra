from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers.kami_assembler import KamiAssembler


def test_create_kami_assembler():
    ka = KamiAssembler()

