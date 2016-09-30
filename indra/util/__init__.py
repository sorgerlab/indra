from __future__ import absolute_import, division, print_function, \
                       unicode_literals
import sys
import csv
import xml.etree.ElementTree as ET

if sys.version_info[0] >= 3:
    non_unicode = bytes
else:
    non_unicode = str

def unicode_strs(obj):
    if isinstance(obj, non_unicode):
        return False
    # Check for an iterable
    if isinstance(obj, list) or isinstance(obj, tuple) or \
       isinstance(obj, set):
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


def decode_obj(obj, encoding='utf-8'):
    if isinstance(obj, non_unicode):
        return obj.decode(encoding)
    elif isinstance(obj, list) or isinstance(obj, tuple):
        return [decode_obj(item) for item in obj]
    elif hasattr(obj, '__dict__'):
        for k, v in obj.__dict__.items():
            obj.__dict__[k] = decode_obj(v)
        return obj
    elif isinstance(obj, dict):
        dec_obj = {}
        for k, v in obj.items():
            dec_k = decode_obj(k)
            dec_v = decode_obj(v)
            dec_obj[dec_k] = dec_v
        return dec_obj
    else:
        return obj


def read_unicode_csv(filename, delimiter=',', quotechar='"',
                     quoting=csv.QUOTE_MINIMAL, lineterminator='\n',
                     encoding='utf-8'):
    # Python 3 version
    if sys.version_info[0] >= 3:
        # Open the file in text mode with given encoding
        # Set newline arg to '' (see https://docs.python.org/3/library/csv.html)
        with open(filename, 'r', newline='', encoding=encoding) as f:
            # Next, get the csv reader, with unicode delimiter and quotechar
            csv_reader = csv.reader(f, delimiter=delimiter, quotechar=quotechar,
                                 quoting=quoting, lineterminator=lineterminator)
            # Now, return the (already decoded) unicode csv_reader generator
            for row in csv_reader:
                yield row
    # Python 2 version
    else:
        # Open the file, no encoding specified
        with open(filename, 'rb') as f:
            # Next, get the csv reader, passing delimiter and quotechar as
            # bytestrings rather than unicode
            csv_reader = csv.reader(f, delimiter=delimiter.encode(encoding),
                                 quotechar=quotechar.encode(encoding),
                                 quoting=quoting, lineterminator=lineterminator)
            # Iterate over the file and decode each string into unicode
            for row in csv_reader:
                yield [cell.decode(encoding) for cell in row]


def write_unicode_csv(filename, rows, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL, lineterminator='\n',
                       encoding='utf-8'):
    # Python 3 version
    if sys.version_info[0] >= 3:
        # Open the file in text mode with given encoding
        # Set newline arg to '' (see https://docs.python.org/3/library/csv.html)
        with open(filename, 'w', newline='', encoding=encoding) as f:
            # Next, get the csv writer, with unicode delimiter and quotechar
            csv_writer = csv.writer(f, delimiter=delimiter, quotechar=quotechar,
                                quoting=quoting, lineterminator=lineterminator)
            # Write the rows to the file
            csv_writer.writerows(rows)
    # Python 2 version
    else:
        # Open the file, no encoding specified
        with open(filename, 'w') as f:
            # Next, get the csv writer, passing delimiter and quotechar as
            # bytestrings rather than unicode
            csv_writer = csv.writer(f, delimiter=delimiter.encode(encoding),
                                quotechar=quotechar.encode(encoding),
                                quoting=quoting, lineterminator=lineterminator)
            for row in rows:
                csv_writer.writerow([cell.encode(encoding) for cell in row])


if sys.version_info[0] >= 3:
    def UnicodeXMLTreeBuilder():
        return None
else:
    class UnicodeXMLTreeBuilder(ET.XMLTreeBuilder):
        # See this thread:
        # http://www.gossamer-threads.com/lists/python/python/728903
        def _fixtext(self, text):
            return text
