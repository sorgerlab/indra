from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
from copy import copy


__all__ = ['Evidence']


import sys
import json
import textwrap
from collections import OrderedDict as _o
from .util import *
from .context import Context


@python_2_unicode_compatible
class Evidence(object):
    """Container for evidence supporting a given statement.

    Parameters
    ----------
    source_api : str or None
        String identifying the INDRA API used to capture the statement,
        e.g., 'trips', 'biopax', 'bel'.
    source_id : str or None
        For statements drawn from databases, ID of the database entity
        corresponding to the statement.
    pmid : str or None
        String indicating the Pubmed ID of the source of the statement.
    text : str
        Natural language text supporting the statement.
    annotations : dict
        Dictionary containing additional information on the
        context of the statement, e.g., species, cell line,
        tissue type, etc. The entries may vary depending on
        the source of the information.
    epistemics : dict
        A dictionary describing various forms of epistemic
        certainty associated with the statement.
    context: Context or None
        A context object
    text_refs : dict
        A dictionary of various reference ids to the source text, e.g.
        DOI, PMID, URL, etc.


    There are some attributes which are not set by the parameters above:

    source_hash : int
        A hash calculated from the evidence text, source api, and pmid and/or
        source_id if available. This is generated automatcially when the object
        is instantiated.
    stmt_tag : int
        This is a hash calculated by a Statement to which this evidence refers,
        and is set by said Statement. It is useful for tracing ownership of
        an Evidence object.
    """
    def __init__(self, source_api=None, source_id=None, pmid=None, text=None,
                 annotations=None, epistemics=None, context=None,
                 text_refs=None):
        self.source_api = source_api
        self.source_id = source_id
        self.pmid = pmid
        self.text_refs = {}
        if pmid is not None:
            self.text_refs['PMID'] = pmid
        if text_refs is not None:
            self.text_refs.update(text_refs)
        self.text = text
        if annotations:
            self.annotations = annotations
        else:
            self.annotations = {}
        if epistemics:
            self.epistemics = epistemics
        else:
            self.epistemics = {}
        self.context = context
        self.source_hash = None
        self.get_source_hash()
        self.stmt_tag = None

    def __setstate__(self, state):
        if 'context' not in state:
            state['context'] = None
        if 'text_refs' not in state:
            state['text_refs'] = {}
        if 'stmt_tag' not in state:
            state['stmt_tag'] = None
        if 'source_hash' not in state:
            state['source_hash'] = None
        self.__dict__ = state

    def get_source_hash(self, refresh=False):
        """Get a hash based off of the source of this statement.

        The resulting value is stored in the source_hash attribute of the class
        and is preserved in the json dictionary.
        """
        if hasattr(self, 'source_hash') and self.source_hash is not None \
                and not refresh:
            return self.source_hash
        s = str(self.source_api) + str(self.source_id)
        if self.text and isinstance(self.text, str):
            s += self.text
        elif self.pmid and isinstance(self.pmid, str):
            s += self.pmid
        self.source_hash = make_hash(s, 16)
        return self.source_hash

    def matches_key(self):
        key_lst = [self.source_api, self.source_id, self.pmid,
                   self.text]
        for d in [self.annotations, self.epistemics]:
            d_key = list(d.items())
            d_key.sort()
            key_lst.append(d_key)
        key = str(key_lst)
        return key.replace('"', '').replace('\'', '').replace('None', '~')[1:-1]

    def equals(self, other):
        matches = (self.source_api == other.source_api) and \
                  (self.source_id == other.source_id) and \
                  (self.pmid == other.pmid) and \
                  (self.text == other.text) and \
                  (self.annotations == other.annotations) and \
                  (self.epistemics == other.epistemics) and \
                  (self.context == other.context)
        return matches

    def to_json(self):
        """Convert the evidence object into a JSON dict."""
        json_dict = _o({})
        if self.source_api:
            json_dict['source_api'] = self.source_api
        if self.pmid:
            json_dict['pmid'] = self.pmid
        if self.source_id:
            json_dict['source_id'] = self.source_id
        if self.text:
            json_dict['text'] = self.text
        if self.annotations:
            json_dict['annotations'] = self.annotations
        if self.epistemics:
            json_dict['epistemics'] = self.epistemics
        if self.context:
            json_dict['context'] = self.context.to_json()
        if self.text_refs:
            json_dict['text_refs'] = self.text_refs
        json_dict['source_hash'] = self.get_source_hash()
        if self.stmt_tag:
            json_dict['stmt_tag'] = self.stmt_tag
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        source_api = json_dict.get('source_api')
        source_id = json_dict.get('source_id')
        pmid = json_dict.get('pmid')
        text = json_dict.get('text')
        annotations = copy(json_dict.get('annotations', {}))
        epistemics = copy(json_dict.get('epistemics', {}))
        context_entry = json_dict.get('context')
        text_refs = copy(json_dict.get('text_refs', {}))
        if context_entry:
            context = Context.from_json(context_entry)
        else:
            context = None
        stmt_tag = json_dict.get('stmt_tag')
        # Note that the source hash will be re-generated upon loading, so if
        # any of the relevant attributes used to create the hash changed, the
        # hash will also have changed.
        ev = Evidence(source_api=source_api, source_id=source_id,
                      pmid=pmid, text=text, annotations=annotations,
                      epistemics=epistemics, context=context,
                      text_refs=text_refs)
        ev.stmt_tag = stmt_tag
        return ev

    def __str__(self):
        ev_str = 'Evidence('
        tab_len = len(ev_str)

        def _indented_join(s_list, depth):
            return '\n'.join(' '*depth + s for s in s_list).lstrip(' ')

        lines = []

        def _add_line(name, s):
            lines.append('%s=%s' % (name, s))

        def _format_line(name, s):
            return _add_line(name, "'%s'" % s)

        def _format_dict(d, name, indent=9):
            s = json.dumps(d, indent=1)
            s = _indented_join(s.splitlines(), indent+len(name)+1)
            return _add_line(name, s)

        if self.source_api:
            _format_line('source_api', self.source_api)
        if self.pmid:
            _format_line('pmid', self.pmid)
        if self.source_id:
            _format_line('source_id', self.source_id)
        if self.text:
            txt = _indented_join(textwrap.wrap(self.text, width=65),
                                 tab_len+6)
            _format_line('text', txt)
        if self.annotations:
            _format_dict(self.annotations, 'annotations')
        if self.context:
            _format_dict(self.context.to_json(), 'context')
        if self.epistemics:
            _format_dict(self.epistemics, 'epistemics')
        if self.text_refs:
            _format_dict(self.text_refs, 'text_refs')

        div = ',\n' + ' '*9
        ev_str += div.join(lines)
        if len(ev_str.splitlines()) > 1:
            ev_str += '\n' + ' '*9
        ev_str += ')\n\n'
        return ev_str

    def __repr__(self):
        if sys.version_info[0] >= 3:
            return str(self)
        else:
            return str(self).encode('utf-8')


