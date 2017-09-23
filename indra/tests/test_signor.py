from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from indra.sources.signor import SignorProcessor

signor_test_path = join(dirname(__file__), '..', '..', 'data',
                        'all_data_23_09_17.csv')

def test_parse_csv():
    sp = SignorProcessor(signor_test_path)

if __name__ == '__main__':
    test_parse_csv()
