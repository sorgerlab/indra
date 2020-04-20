__all__ = ['get_valid_residue', 'get_valid_location', 'activity_types',
           'amino_acids', 'amino_acids_reverse',
           'InvalidLocationError', 'InvalidResidueError']


import os
import rdflib


def get_valid_residue(residue):
    """Check if the given string represents a valid amino acid residue."""
    if residue is not None and amino_acids.get(residue) is None:
        res = amino_acids_reverse.get(residue.lower())
        if res is None:
            raise InvalidResidueError(residue)
        else:
            return res
    return residue


def get_valid_location(location):
    """Check if the given location represents a valid cellular component."""
    from indra.databases import go_client
    # If we're given None, return None
    if location is None:
        return None
    # Otherwise, er look up a GO ID by name or synonym and then get the
    # canonical name for that ID. If the location provided is not a valid
    # label of synonym, we raise an error
    go_id = go_client.get_go_id_from_label_or_synonym(location)
    if not go_id:
        raise InvalidLocationError(location)
    standard_name = go_client.get_go_label(go_id)
    return standard_name


def _read_activity_types():
    """Read types of valid activities from a resource file."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ac_file = os.path.join(this_dir, os.pardir, 'resources',
                           'activity_hierarchy.rdf')
    g = rdflib.Graph()
    with open(ac_file, 'r'):
        g.parse(ac_file, format='nt')
    act_types = set()
    for s, _, o in g:
        subj = s.rpartition('/')[-1]
        obj = o.rpartition('/')[-1]
        act_types.add(subj)
        act_types.add(obj)
    return sorted(list(act_types))


activity_types = _read_activity_types()


def _read_amino_acids():
    """Read the amino acid information from a resource file."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    aa_file = os.path.join(this_dir, os.pardir, 'resources', 'amino_acids.tsv')
    amino_acids = {}
    amino_acids_reverse = {}
    with open(aa_file, 'rt') as fh:
        lines = fh.readlines()
    for lin in lines[1:]:
        terms = lin.strip().split('\t')
        key = terms[2]
        val = {'full_name': terms[0],
               'short_name': terms[1],
               'indra_name': terms[3]}
        amino_acids[key] = val
        for v in val.values():
            amino_acids_reverse[v] = key
    return amino_acids, amino_acids_reverse


amino_acids, amino_acids_reverse = _read_amino_acids()


class InvalidResidueError(ValueError):
    """Invalid residue (amino acid) name."""
    def __init__(self, name):
        ValueError.__init__(self, "Invalid residue name: '%s'" % name)


class InvalidLocationError(ValueError):
    """Invalid cellular component name."""
    def __init__(self, name):
        ValueError.__init__(self, "Invalid location name: '%s'" % name)
