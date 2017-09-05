"""
Import this file in your script if you want to ensure that it
is only used in python 3.
"""
import sys
if sys.version_info.major < 3:
    raise Exception('This should be used in python 3 only.')
