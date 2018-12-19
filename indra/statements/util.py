from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible


__all__ = ['make_hash']


from hashlib import md5


def make_hash(s, n_bytes):
    """Make the hash from a matches key."""
    raw_h = int(md5(s.encode('utf-8')).hexdigest()[:n_bytes], 16)
    # Make it a signed int.
    return 16**n_bytes//2 - raw_h