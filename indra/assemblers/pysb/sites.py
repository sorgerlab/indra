__all__ = ['abbrevs', 'states', 'mod_acttype_map', 'get_binding_site_name',
           'get_mod_site_name']
from indra.statements import *
from indra.ontology.bio import bio_ontology
from .common import _n

abbrevs = {
    'phosphorylation': 'phospho',
    'ubiquitination': 'ub',
    'farnesylation': 'farnesyl',
    'hydroxylation': 'hydroxyl',
    'acetylation': 'acetyl',
    'sumoylation': 'sumo',
    'glycosylation': 'glycosyl',
    'methylation': 'methyl',
    'ribosylation': 'ribosyl',
    'geranylgeranylation': 'geranylgeranyl',
    'palmitoylation': 'palmitoyl',
    'myristoylation': 'myryl',
    'modification': 'mod',
}

states = {
    'phosphorylation': ['u', 'p'],
    'ubiquitination': ['n', 'y'],
    'farnesylation': ['n', 'y'],
    'hydroxylation': ['n', 'y'],
    'acetylation': ['n', 'y'],
    'sumoylation': ['n', 'y'],
    'glycosylation': ['n', 'y'],
    'methylation': ['n', 'y'],
    'geranylgeranylation': ['n', 'y'],
    'palmitoylation': ['n', 'y'],
    'myristoylation': ['n', 'y'],
    'ribosylation': ['n', 'y'],
    'modification': ['n', 'y'],
}

mod_acttype_map = {
    Phosphorylation: 'kinase',
    Dephosphorylation: 'phosphatase',
    Hydroxylation: 'catalytic',
    Dehydroxylation: 'catalytic',
    Sumoylation: 'catalytic',
    Desumoylation: 'catalytic',
    Acetylation: 'catalytic',
    Deacetylation: 'catalytic',
    Glycosylation: 'catalytic',
    Deglycosylation: 'catalytic',
    Ribosylation: 'catalytic',
    Deribosylation: 'catalytic',
    Ubiquitination: 'catalytic',
    Deubiquitination: 'catalytic',
    Farnesylation: 'catalytic',
    Defarnesylation: 'catalytic',
    Palmitoylation: 'catalytic',
    Depalmitoylation: 'catalytic',
    Myristoylation: 'catalytic',
    Demyristoylation: 'catalytic',
    Geranylgeranylation: 'catalytic',
    Degeranylgeranylation: 'catalytic',
    Methylation: 'catalytic',
    Demethylation: 'catalytic',
}


def get_binding_site_name(agent):
    """Return a binding site name from a given agent."""
    # Try to construct a binding site name based on parent
    grounding = agent.get_grounding()
    if grounding != (None, None):
        top_parents = bio_ontology.get_top_level_parents(*grounding)
        if top_parents:
            parent_name = bio_ontology.get_name(*top_parents[0])
            if parent_name:
                return _n(parent_name).lower()
    return _n(agent.name).lower()


def get_mod_site_name(mod_condition):
    """Return site names for a modification."""
    if mod_condition.residue is None:
        mod_str = abbrevs[mod_condition.mod_type]
    else:
        mod_str = mod_condition.residue
    mod_pos = mod_condition.position if \
        mod_condition.position is not None else ''
    name = ('%s%s' % (mod_str, mod_pos))
    return name
