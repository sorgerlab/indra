__all__ = ['_run_eidos_on_text', 'process_text_bio',
           'process_json_bio', 'process_json_bio_entities',
           'process_text_bio_entities', 'eidos_reader',
           'initialize_reader']

import json
import logging
from indra.sources.eidos import client as eidos_client
from .bio_processor import EidosBioProcessor

logger = logging.getLogger(__name__)


try:
    # For text reading
    from .reader import EidosReader
    eidos_reader = EidosReader()
except Exception as e:
    logger.warning('Could not instantiate Eidos reader, local reading '
                   'will not be available.')
    eidos_reader = None


def _run_eidos_on_text(text, save_json='eidos_output.json',
                       webservice=None):
    if not webservice:
        if eidos_reader is None:
            logger.error('Eidos reader is not available.')
            return None
        json_dict = eidos_reader.process_text(text)
    else:
        if webservice.endswith('/'):
            webservice = webservice[:-1]
        json_dict = eidos_client.process_text(text, webservice=webservice)
    if json_dict and save_json:
        with open(save_json, 'wt') as fh:
            json.dump(json_dict, fh, indent=2)
    return json_dict


def process_text_bio(text, save_json='eidos_output.json', webservice=None,
                     grounder=None):
    """Return an EidosProcessor by processing the given text.

    This constructs a reader object via Java and extracts mentions
    from the text. It then serializes the mentions into JSON and
    processes the result with process_json.

    Parameters
    ----------
    text : str
        The text to be processed.
    save_json : Optional[str]
        The name of a file in which to dump the JSON output of Eidos.
    webservice : Optional[str]
        An Eidos reader web service URL to send the request to.
        If None, the reading is assumed to be done with the Eidos JAR rather
        than via a web service. Default: None
    grounder : Optional[function]
        A function which takes a text and an optional context as argument
        and returns a dict of groundings.

    Returns
    -------
    ep : EidosProcessor
        An EidosProcessor containing the extracted INDRA Statements in its
        statements attribute.
    """
    json_dict = _run_eidos_on_text(text, save_json, webservice)
    if json_dict:
        return process_json_bio(json_dict, grounder=grounder)
    return None


def process_json_bio(json_dict, grounder=None):
    """Return EidosProcessor with grounded Activation/Inhibition statements.

    Parameters
    ----------
    json_dict : dict
        The JSON-LD dict to be processed.
    grounder : Optional[function]
        A function which takes a text and an optional context as argument
        and returns a dict of groundings.

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in its statements attribute.
    """
    from indra.sources.eidos.bio_processor import EidosBioProcessor
    ep = EidosBioProcessor(json_dict, grounder=grounder)
    ep.extract_statements()
    return ep


def process_json_bio_entities(json_dict, grounder=None, with_coords=False):
    """Return INDRA Agents grounded to biological ontologies extracted
    from Eidos JSON-LD.

    Parameters
    ----------
    json_dict : dict
        The JSON-LD dict to be processed.
    grounder : Optional[function]
        A function which takes a text and an optional context as argument
        and returns a dict of groundings.
    with_coords : Optional[bool]
        If True, the Agents will have their coordinates returned along
        with them in a tuple. Default: False

    Returns
    -------
    list of indra.statements.Agent
        A list of INDRA Agents which are derived from concepts extracted
        by Eidos from text.
    """
    from .bio_processor import get_agent_bio
    if not json_dict:
        return []
    ep = EidosBioProcessor(json_dict, grounder=grounder)
    ep.extract_causal_relations()
    ep.extract_events()
    events = ep.get_all_events()
    agents = []
    for event in events:
        context = event.evidence[0].text
        agent = get_agent_bio(event.concept, context=context,
                              grounder=grounder)
        if with_coords:
            prov = event.evidence[0].annotations['provenance']
            pos = prov[0]['documentCharPositions']
            start_coord, end_coord = pos[0]['start'], pos[0]['end']
            # We make sure this is a proper coordinate, in which case
            # we add one to ensure consistency with coordinate boundaries
            # from other INDRA sources
            if isinstance(end_coord, int):
                end_coord += 1
            agents.append((agent, (start_coord, end_coord)))
        else:
            agents.append(agent)
    return agents


def process_text_bio_entities(text, webservice=None, grounder=None):
    """Return INDRA Agents grounded to biological ontologies extracted
    from text.

    Parameters
    ----------
    text : str
        Text to be processed.
    webservice : Optional[str]
        An Eidos reader web service URL to send the request to.
        If None, the reading is assumed to be done with the Eidos JAR rather
        than via a web service. Default: None
    grounder : Optional[function]
        A function which takes a text and an optional context as argument
        and returns a dict of groundings.

    Returns
    -------
    list of indra.statements.Agent
        A list of INDRA Agents which are derived from concepts extracted
        by Eidos from text.
    """
    json_dict = _run_eidos_on_text(text, None, webservice=webservice)
    return process_json_bio_entities(json_dict, grounder=grounder)


def initialize_reader():
    """Instantiate an Eidos reader for fast subsequent reading."""
    eidos_reader.process_text('')
