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
