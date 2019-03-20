__all__ = ['perm_cache']

import json
import pickle
from functools import update_wrapper
from os.path import exists


def perm_cache(cache_type='pkl', cache_file=None):
    class PermCache(object):
        _cache_type = cache_type
        _cache_file = cache_file

        def __init__(self, func):
            if self._cache_type not in ['pkl', 'json']:
                raise ValueError("Invalid cache type: %s" % self._cache_type)
            self._cache_type = self._cache_type
            self.func = func
            if self._cache_file is None:
                self._cache_file = (func.__code__.co_filename
                                    .replace('.py', '.' + self.func.__name__))
                self._cache_file += '.cache'
            if self._cache_file.endswith('.py'):
                self._cache_file = self._cache_file.replace('.py',
                                                          '.' + self._cache_type)
            else:
                self._cache_file += '.' + self._cache_type

            if exists(self._cache_file):
                if self._cache_type == 'pkl':
                    with open(self._cache_file, 'rb') as f:
                        self.cache = pickle.load(f)
                elif self._cache_type == 'json':
                    with open(self._cache_file, 'r') as f:
                        self.cache = json.load(f)
            else:
                self.cache = {}
            self.__cache_info = dict.fromkeys(['added', 'read', 'total'], 0)
            update_wrapper(self, func)
            return

        def __call__(self, *args, **kwargs):
            key = ' '.join(args) \
                  + ' '.join(['%s=%s' % (k, v) for k, v in kwargs.items()])
            self.__cache_info['total'] += 1
            try:
                res = self.cache[key]
                self.__cache_info['read'] += 1
            except KeyError:
                res = self.func(*args, **kwargs)
                self.cache[key] = res
                self.__cache_info['added'] += 1
            return res

        def cache_info(self):
            return self.__cache_info.copy()

        def stash_cache(self):
            if self._cache_type == 'pkl':
                with open(self._cache_file, 'wb') as f:
                    pickle.dump(self.cache, f)
            elif self._cache_type == 'json':
                with open(self._cache_file, 'w') as f:
                    json.dump(self.cache, f, indent=2)
            return

    return PermCache
