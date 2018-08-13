from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str


class NestedDict(dict):
    """A dict-like object that recursively populates elements of a dict."""

    def __getitem__(self, key):
        if key not in self.keys():
            val = self.__class__()
            self.__setitem__(key, val)
        else:
            val = dict.__getitem__(self, key)
        return val

    def __repr__(self):
        sub_str = dict.__repr__(self)[1:-1]
        if not sub_str:
            return self.__class__.__name__ + '()'
        # This does not completely generalize, but it works for most cases.
        for old, new in [('), ', '),\n'), ('\n', '\n  ')]:
            sub_str = sub_str.replace(old, new)
        return'%s(\n  %s\n)' % (self.__class__.__name__, sub_str)

    def __str__(self):
        return self.__repr__()

    def export_dict(self):
        "Convert this into an ordinary dict (of dicts)."
        return {k: v.export_dict() if isinstance(v, self.__class__) else v
                for k, v in self.items()}

    def get(self, key):
        "Find the first value within the tree which has the key."
        if key in self.keys():
            return self[key]
        else:
            res = None
            for v in self.values():
                # This could get weird if the actual expected returned value
                # is None, especially in teh case of overlap. Any ambiguity
                # would be resolved by get_path(s).
                if hasattr(v, 'get'):
                    res = v.get(key)
                if res is not None:
                    break
            return res

    def get_path(self, key):
        "Like `get`, but also return the path taken to the value."
        if key in self.keys():
            return (key,), self[key]
        else:
            key_path, res = (None, None)
            for sub_key, v in self.items():
                if isinstance(v, self.__class__):
                    key_path, res = v.get_path(key)
                elif hasattr(v, 'get'):
                    res = v.get(key)
                    key_path = (key,) if res is not None else None
                if res is not None and key_path is not None:
                    key_path = (sub_key,) + key_path
                    break
            return key_path, res

    def gets(self, key):
        "Like `get`, but return all matches, not just the first."
        result_list = []
        if key in self.keys():
            result_list.append(self[key])
        for v in self.values():
            if isinstance(v, self.__class__):
                sub_res_list = v.gets(key)
                for res in sub_res_list:
                    result_list.append(res)
            elif isinstance(v, dict):
                if key in v.keys():
                    result_list.append(v[key])
        return result_list

    def get_paths(self, key):
        "Like `gets`, but include the paths, like `get_path` for all matches."
        result_list = []
        if key in self.keys():
            result_list.append(((key,), self[key]))
        for sub_key, v in self.items():
            if isinstance(v, self.__class__):
                sub_res_list = v.get_paths(key)
                for key_path, res in sub_res_list:
                    result_list.append(((sub_key,) + key_path, res))
            elif isinstance(v, dict):
                if key in v.keys():
                    result_list.append(((sub_key, key), v[key]))
        return result_list

    def get_leaves(self):
        """Get the deepest entries as a flat set."""
        ret_set = set()
        for val in self.values():
            if isinstance(val, self.__class__):
                ret_set |= val.get_leaves()
            elif isinstance(val, dict):
                ret_set |= set(val.values())
            elif isinstance(val, list):
                ret_set |= set(val)
            elif isinstance(val, set):
                ret_set |= val
            else:
                ret_set.add(val)
        return ret_set
