from sys import version_info
from functools import wraps
from nose import SkipTest

IS_PY3 = True
if version_info.major is not 3:
    IS_PY3 = False


def needs_py3(func):
    @wraps(func)
    def test_with_py3_func(*args, **kwargs):
        if not IS_PY3:
            raise SkipTest("This tests features only supported in Python 3.x")
        return func(*args, **kwargs)
    return test_with_py3_func


def skip_if(condition, reason=None):
    """Skip test if condition is true"""
    def decorate(func):
        wraps(func)
        func_name = func.__name__

        def f(*args, **kwargs):
            if not condition:
                raise SkipTest("'%s' skipped: %s" % (func_name, reason))
            else:
                return func(*args, **kwargs)
        return f
    return decorate
