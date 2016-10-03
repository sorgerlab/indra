from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from rasmachine import make_status_message
from rasmachine import _increment_ndex_ver

stats = {}
stats['new_abstracts'] = 0
stats['new_papers'] = 0
stats['orig_stmts'] = 10
stats['new_stmts'] = 10
stats['orig_top'] = 10
stats['new_top'] = 10

def test_ndex_ver():
    assert(_increment_ndex_ver(None) == '1.0')
    assert(_increment_ndex_ver('') == '1.0')
    assert(_increment_ndex_ver('1.0') == '1.1')
    assert(_increment_ndex_ver('1.9') == '1.10')
    assert(_increment_ndex_ver('2.10') == '2.11')

def test_noabs_nopaper():
    s = stats.copy()
    status_msg = make_status_message(s)
    assert(status_msg is None)

def test_absonly():
    s = stats.copy()
    s['new_abstracts'] = 7
    s['new_top'] = 12
    status_msg = make_status_message(s)
    assert(status_msg == 'Today I read 7 abstracts, ' +\
                         'and learned 2 new mechanisms!')

def test_papersonly():
    s = stats.copy()
    s['new_papers'] = 3
    s['new_top'] = 12
    status_msg = make_status_message(s)
    assert(status_msg == 'Today I read 3 papers, ' +\
                         'and learned 2 new mechanisms!')

def test_abs_and_paper():
    s = stats.copy()
    s['new_papers'] = 1
    s['new_abstracts'] = 1
    s['new_top'] = 12
    status_msg = make_status_message(s)
    assert(status_msg == 'Today I read 1 paper and 1 abstract, ' +\
                         'and learned 2 new mechanisms!')

def test_abs_and_papers():
    s = stats.copy()
    s['new_papers'] = 2
    s['new_abstracts'] = 2
    s['new_top'] = 12
    status_msg = make_status_message(s)
    assert(status_msg == 'Today I read 2 papers and 2 abstracts, ' +\
                         'and learned 2 new mechanisms!')
