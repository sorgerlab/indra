from __future__ import absolute_import, division, print_function, \
                       unicode_literals
import sys

if sys.version_info > (3, 0):
    non_unicode = bytes
else:
    non_unicode = str

def unicode_strs(obj):
    if isinstance(obj, non_unicode):
        return False
    # Check for an iterable
    if isinstance(obj, list) or isinstance(obj, tuple):
        for item in obj:
            has_unicode_strs = unicode_strs(item)
            if not has_unicode_strs:
                return False
    if hasattr(obj, '__dict__'):
        for item in obj.__dict__.values():
            has_unicode_strs = unicode_strs(item)
            if not has_unicode_strs:
                return False
    if isinstance(obj, dict):
        for k, v in obj.items():
            k_has_unicode_strs = unicode_strs(k)
            v_has_unicode_strs = unicode_strs(v)
            if not k_has_unicode_strs or not v_has_unicode_strs:
                return False
    return True


#foo = {u'a':u'b', u'c':[u'a', u'b', u'c']}
#foo = {'a':'b', 'c':['a', 'b', 'c']}
#print(unicode_strs(foo))
