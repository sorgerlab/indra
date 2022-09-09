from collections import defaultdict

from indra.util import read_unicode_csv
from indra.resources import get_resource_path


def get_id_from_name(name):
    return mgi_name_to_id.get(name)


def get_name_from_id(mgi_id):
    return mgi_id_to_name.get(mgi_id)


def get_synonyms(mgi_id):
    return mgi_synonyms.get(mgi_id)


def get_id_from_name_synonym(name_synonym):
    mgi_id = mgi_name_to_id.get(name_synonym)
    if mgi_id:
        return mgi_id
    mgi_ids = mgi_synonyms_reverse.get(name_synonym)
    if mgi_ids:
        if len(mgi_ids) == 1:
            return mgi_ids[0]
        else:
            return mgi_ids
    return None


def _read_mgi():
    fname = get_resource_path('mgi_entries.tsv')
    mgi_id_to_name = {}
    mgi_name_to_id = {}
    mgi_synonyms = {}
    mgi_synonyms_reverse = defaultdict(list)
    for mgi_id, name, synonyms_str in read_unicode_csv(fname, '\t'):
        if name:
            mgi_id_to_name[mgi_id] = name
            mgi_name_to_id[name] = mgi_id
        if synonyms_str:
            synonyms = synonyms_str.split('|')
            mgi_synonyms[mgi_id] = synonyms
            for synonym in synonyms:
                mgi_synonyms_reverse[synonym] = mgi_id

    return mgi_id_to_name, mgi_name_to_id, mgi_synonyms, \
        dict(mgi_synonyms_reverse)


mgi_id_to_name, mgi_name_to_id, mgi_synonyms, mgi_synonyms_reverse = _read_mgi()
