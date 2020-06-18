"""Copied from here
https://bugs.python.org/issue13831#msg187247

Only change to code is formatting and use of f-strings
"""
import sys
import traceback


class WrapException(Exception):
    """Exception wrapper for multiprocessing"""
    def __init__(self):
        exc_type, exc_value, exc_tb = sys.exc_info()
        self.exception = exc_value
        self.formatted = ''.join(traceback.format_exception(exc_type,
                                                            exc_value,
                                                            exc_tb))

    def __str__(self):
        return f'{Exception.__str__(self)}\n' \
               f'Original traceback:\n{self.formatted}'
