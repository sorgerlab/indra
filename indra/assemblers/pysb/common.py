__all__ = ['_n']
import re


def _n(name):
    """Return valid PySB name."""
    n = name.encode('ascii', errors='ignore').decode('ascii')
    n = re.sub('[^A-Za-z0-9_]', '_', n)
    n = re.sub(r'(^[0-9].*)', r'p\1', n)
    return n
