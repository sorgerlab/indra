__all__ = ['get_valid_residue', 'activity_types',
           'amino_acids', 'amino_acids_reverse',
           'InvalidLocationError', 'InvalidResidueError']

import os


def get_valid_residue(residue):
    """Check if the given string represents a valid amino acid residue."""
    if residue is not None and amino_acids.get(residue) is None:
        res = amino_acids_reverse.get(residue.lower())
        if res is None:
            raise InvalidResidueError(residue)
        else:
            return res
    return residue


# FIXME: this list could be read from the ontology or another resource file
activity_types = {'transcription', 'activity', 'catalytic' 'gtpbound',
                  'kinase', 'phosphatase', 'gef', 'gap'}


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
