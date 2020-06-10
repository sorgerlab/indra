"""Module containing the implementation of an IndraOntology for the
 general biology use case."""
__all__ = ['bio_ontology', 'BioOntology']

from indra.config import get_config
from .ontology import BioOntology
from ..virtual import VirtualOntology

indra_ontology_url = get_config('INDRA_ONTOLOGY_URL')
bio_ontology = BioOntology() if not indra_ontology_url else \
    VirtualOntology(url=indra_ontology_url)