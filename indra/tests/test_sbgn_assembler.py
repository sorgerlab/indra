from os.path import dirname, join
from hashlib import sha1
from indra.assemblers.sbgn_assembler import text_to_sbgn

test_small_file = join(dirname(__file__), 'test_small.xml')

# Simple smoke test to ensure the test_small file produces known-correct output.
def test_assembler():
    sbgn_output = text_to_sbgn(trips_xml=open(test_small_file).read())
    assert (sha1(sbgn_output).hexdigest() ==
            'cf5436f3687db9ff61da73af76595b3836dd65de')
