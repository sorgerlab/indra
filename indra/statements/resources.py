from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str


__all__ = ['get_valid_residue', 'get_valid_location', 'activity_types',
           'cellular_components', 'cellular_components_reverse',
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
    # If we're given None, return None
    if location is not None and cellular_components.get(location) is None:
        loc = cellular_components_reverse.get(location)
        if loc is None:
            raise InvalidLocationError(location)
        else:
            return loc
    return location


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


def _read_cellular_components():
    """Read cellular components from a resource file."""
    # Here we load a patch file in addition to the current cellular components
    # file to make sure we don't error with InvalidLocationError with some
    # deprecated cellular location names
    this_dir = os.path.dirname(os.path.abspath(__file__))
    cc_file = os.path.join(this_dir, os.pardir, 'resources',
                           'cellular_components.tsv')
    cc_patch_file = os.path.join(this_dir, os.pardir, 'resources',
                                 'cellular_components_patch.tsv')
    cellular_components = {}
    cellular_components_reverse = {}
    with open(cc_file, 'rt') as fh:
        lines = list(fh.readlines())
    # We add the patch to the end of the lines list
    with open(cc_patch_file, 'rt') as fh:
        lines += list(fh.readlines())
    for lin in lines[1:]:
        terms = lin.strip().split('\t')
        cellular_components[terms[1]] = terms[0]
        # If the GO -> name mapping doesn't exist yet, we add a mapping
        # but if it already exists (i.e. the try doesn't error) then
        # we don't add the GO -> name mapping. This ensures that names from
        # the patch file aren't mapped to in the reverse list.
        try:
            cellular_components_reverse[terms[0]]
        except KeyError:
            cellular_components_reverse[terms[0]] = terms[1]
    return cellular_components, cellular_components_reverse


cellular_components, cellular_components_reverse = _read_cellular_components()


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
