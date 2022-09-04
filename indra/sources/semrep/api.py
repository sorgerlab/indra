__all__ = ['process_xml_file']

from xml.etree import ElementTree as ET
from .processor import SemRepXmlProcessor


def process_xml_file(fname, use_gilda_grounding=False, predicate_mappings=None):
    """Process a SemRep output XML file and extract INDRA Statements.

    Parameters
    ----------
    fname : str
        The name of the SemRep output XML file.
    use_gilda_grounding : Optional[bool]
        If True, Gilda is used to re-ground entities and assing identifiers.
        Default: False
    predicate_mappings : Optional[dict]
        Allows providing a custom mapping of SemRep predicates to INDRA
        Statement types. If not provided, default ones are used.

    Returns
    -------
    SemRepXmlProcessor
        An instance of a SemRepXmlProcessor that carries extracted INDRA
        Statements in its statements attribute.
    """
    tree = ET.parse(fname)
    sp = SemRepXmlProcessor(tree, use_gilda_grounding=use_gilda_grounding,
                            predicate_mappings=predicate_mappings)
    sp.process_statements()
    return sp
