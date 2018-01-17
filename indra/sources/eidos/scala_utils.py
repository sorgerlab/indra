from jnius import autoclass
import logging

logger = logging.getLogger('eidos')

def get_python_list(scala_list):
    """Return list from elements of scala.collection.immutable.List"""
    python_list = []
    for i in range(scala_list.length()):
        python_list.append(scala_list.apply(i))
    return python_list


def get_python_dict(scala_map):
    """Return a dict from entries in a scala.collection.immutable.Map"""
    python_dict = {}
    keys = get_python_list(scala_map.keys().toList())
    for key in keys:
        python_dict[key] = scala_map.apply(key)
    return python_dict


def get_python_json(scala_json):
    """Return a JSON dict from a org.json4s.JsonAST"""
    def convert_node(node):
        if node.__class__.__name__ in ('org.json4s.JsonAST$JValue',
                                       'org.json4s.JsonAST$JObject'):
            # Make a dictionary and then convert each value
            values_raw = get_python_dict(node.values())
            values = {}
            for k, v in values_raw.items():
                values[k] = convert_node(v)
            return values
        elif node.__class__.__name__.startswith('scala.collection.immutable.Map') or \
            node.__class__.__name__ == \
                'scala.collection.immutable.HashMap$HashTrieMap':
            values_raw = get_python_dict(node)
            values = {}
            for k, v in values_raw.items():
                values[k] = convert_node(v)
            return values
        elif node.__class__.__name__ == 'org.json4s.JsonAST$JArray':
            entries_raw = get_python_list(node.values())
            entries = []
            for entry in entries_raw:
                entries.append(convert_node(entry))
            return entries
        elif node.__class__.__name__ == 'scala.collection.immutable.$colon$colon':
            entries_raw = get_python_list(node)
            entries = []
            for entry in entries_raw:
                entries.append(convert_node(entry))
            return entries
        elif node.__class__.__name__ == 'scala.math.BigInt':
            return node.intValue()
        elif node.__class__.__name__ == 'scala.None$':
            return None
        elif node.__class__.__name__ == 'scala.collection.immutable.Nil$':
            return []
        elif isinstance(node, (str, int, float)):
            return node
        else:
            logger.error('Cannot convert %s into Python' %
                         node.__class__.__name__)
            return node.__class__.__name__

    python_json = convert_node(scala_json)
    return python_json

