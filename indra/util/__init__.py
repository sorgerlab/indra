from __future__ import absolute_import, division, print_function, \
                       unicode_literals
import sys
import csv

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

def unicode_csv_reader(f, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL, lineterminator='\n',
                       encoding='utf-8'):
    # Unicode literals will fail in Python 2
    try:
        csv_reader = csv.reader(f, delimiter=delimiter, quotechar=quotechar,
                                quoting=quoting, lineterminator=lineterminator)
    except TypeError:
        csv_reader = csv.reader(f, delimiter=delimiter.encode('utf-8'),
                                quotechar=quotechar.encode('utf-8'),
                                quoting=quoting, lineterminator=lineterminator)
    for row in csv_reader:
        yield [cell.decode(encoding) for cell in row]

def write_unicode_csv(f, rows, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL, lineterminator='\n',
                       encoding='utf-8'):
    try:
        csv_writer = csv.writer(f, delimiter=delimiter, quotechar=quotechar,
                                quoting=quoting, lineterminator=lineterminator)
    except TypeError:
        csv_writer = csv.writer(f, delimiter=delimiter.encode('utf-8'),
                                quotechar=quotechar.encode('utf-8'),
                                quoting=quoting, lineterminator=lineterminator)
    for row in rows:
        csv_writer.writerow([cell.encode(encoding) for cell in row])
    return

#foo = {u'a':u'b', u'c':[u'a', u'b', u'c']}
#foo = {'a':'b', 'c':['a', 'b', 'c']}
#print(unicode_strs(foo))
