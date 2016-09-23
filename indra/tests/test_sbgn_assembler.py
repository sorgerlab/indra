from __future__ import print_function, unicode_literals
from os.path import dirname, join
from hashlib import sha1
from indra.assemblers.sbgn_assembler import text_to_sbgn

test_small_file = join(dirname(__file__), 'test_small.xml')

# Simple smoke test to ensure the test_small file produces known-correct output.
def test_assembler():
    sbgn_output = text_to_sbgn(trips_xml=open(test_small_file).read())
    assert (sha1(sbgn_output).hexdigest() ==
            '99f74446e6e5f0be4a37666c7ac6b74afe5290d0')

