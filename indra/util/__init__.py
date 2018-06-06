from __future__ import absolute_import, division, print_function, \
                       unicode_literals
import sys
import csv
import gzip
import zlib
from io import BytesIO
import xml.etree.ElementTree as ET
try:  # Python 3
    from itertools import zip_longest
except ImportError:  # Python 2
    from itertools import izip_longest as zip_longest


if sys.version_info[0] >= 3:
    non_unicode = bytes
    import pickle
else:
    non_unicode = str
    import cPickle as pickle


def unicode_strs(obj, attr_filter=None):
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
        for item_name, item in obj.__dict__.items():
            if attr_filter and item_name in attr_filter:
                continue
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
                     encoding='utf-8', skiprows=0):
    # Python 3 version
    if sys.version_info[0] >= 3:
        # Open the file in text mode with given encoding
        # Set newline arg to '' (see https://docs.python.org/3/library/csv.html)
        with open(filename, 'r', newline='', encoding=encoding) as f:
            generator = read_unicode_csv_fileobj(f, delimiter=delimiter,
                                            quotechar=quotechar,
                                            quoting=quoting,
                                            lineterminator=lineterminator,
                                            encoding=encoding,
                                            skiprows=skiprows)
            for row in generator:
                yield row
    # Python 2 version
    else:
        # Open the file in binary mode
        with open(filename, 'rb') as f:
            generator = read_unicode_csv_fileobj(f, delimiter=delimiter,
                                            quotechar=quotechar,
                                            quoting=quoting,
                                            lineterminator=lineterminator,
                                            encoding=encoding,
                                            skiprows=skiprows)
            for row in generator:
                yield row


def read_unicode_csv_fileobj(fileobj, delimiter=',', quotechar='"',
                             quoting=csv.QUOTE_MINIMAL, lineterminator='\n',
                             encoding='utf-8', skiprows=0):
    """fileobj can be a StringIO in Py3, but should be a BytesIO in Py2."""
    # Python 3 version
    if sys.version_info[0] >= 3:
        # Next, get the csv reader, with unicode delimiter and quotechar
        csv_reader = csv.reader(fileobj, delimiter=delimiter,
                                quotechar=quotechar, quoting=quoting,
                                lineterminator=lineterminator)
        # Now, return the (already decoded) unicode csv_reader generator
        # Skip rows if necessary
        for skip_ix in range(skiprows):
            next(csv_reader)
        for row in csv_reader:
            yield row
    # Python 2 version
    else:
        # Next, get the csv reader, passing delimiter and quotechar as
        # bytestrings rather than unicode
        csv_reader = csv.reader(fileobj, delimiter=delimiter.encode(encoding),
                             quotechar=quotechar.encode(encoding),
                             quoting=quoting, lineterminator=lineterminator)
        # Iterate over the file and decode each string into unicode
        # Skip rows if necessary
        for skip_ix in range(skiprows):
            next(csv_reader)
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
                csv_writer.writerow([unicode(cell).encode(encoding)
                                     for cell in row])


def zip_string(content, name='gzipped_object'):
    buf = BytesIO()
    gzf = gzip.GzipFile(name, 'wb', 6, buf)
    gzf.write(content.encode('utf8'))
    gzf.close()
    return buf.getvalue()


def unzip_string(gz_obj):
    # Get the content from the object
    gz_body = gz_obj['Body'].read()
    # Decode the gzipped content
    content = zlib.decompress(gz_body, 16+zlib.MAX_WBITS)
    return content.decode('utf8')


if sys.version_info[0] >= 3:
    def UnicodeXMLTreeBuilder():
        return None
else:
    class UnicodeXMLTreeBuilder(ET.XMLTreeBuilder):
        # See this thread:
        # http://www.gossamer-threads.com/lists/python/python/728903
        def _fixtext(self, text):
            return text


def fast_deepcopy(obj):
    """This is a faster implementation of deepcopy via pickle.

    It is meant primarily for sets of Statements with complex hierarchies
    but can be used for any object.
    """
    with BytesIO() as buf:
        pickle.dump(obj, buf)
        buf.seek(0)
        obj_new = pickle.load(buf)
    return obj_new


def lmap(f, xs):
    """A non-lazy version of map."""
    return list(map(f, xs))


def flatten(l):
    """Flatten a nested list."""
    return sum(map(flatten, l), []) \
        if isinstance(l, list) or isinstance(l, tuple) else [l]


def flatMap(f, xs):
    """Map a function onto an iterable and flatten the result."""
    return flatten(lmap(f, xs))


def batch_iter(iterator, batch_size, return_func=None, padding=None):
    """Break an iterable into batches of size batch_size

    Note that `padding` should be set to something (anything) which is NOT a
    valid member of the iterator. For example, None works for [0,1,2,...10], but
    not for ['a', None, 'c', 'd'].

    Parameters
    ----------
    iterator : iterable
        A python object which is iterable.
    batch_size : int
        The size of batches you wish to produce from the iterator.
    return_func : executable or None
        Pass a function that takes a generator and returns an iterable (e.g.
        `list` or `set`). If None, a generator will be returned.
    padding : anything
        This is used internally to ensure that the remainder of the list is
        included. This MUST NOT be a valid element of the iterator.

    Returns
    -------
    An iterator over lists or generators, depending on `return_lists`.
    """
    for batch in zip_longest(*[iter(iterator)]*batch_size, fillvalue=padding):
        gen = (thing for thing in batch if thing is not padding)
        if return_func is None:
            yield gen
        else:
            yield return_func(gen)