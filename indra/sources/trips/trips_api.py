from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from .processor import TripsProcessor
from indra.sources.trips import trips_client

logger = logging.getLogger('trips')

try:
    from .drum_reader import DrumReader
    offline_reading = True
    offline_err = None
except Exception as e:
    offline_err = e
    offline_reading = False


def process_text(text, save_xml_name='trips_output.xml', save_xml_pretty=True,
                 offline=False, service_endpoint='drum'):
    """Return a TripsProcessor by processing text.

    Parameters
    ----------
    text : str
        The text to be processed.
    save_xml_name : Optional[str]
        The name of the file to save the returned TRIPS extraction knowledge
        base XML. Default: trips_output.xml
    save_xml_pretty : Optional[bool]
        If True, the saved XML is pretty-printed. Some third-party tools
        require non-pretty-printed XMLs which can be obtained by setting this
        to False. Default: True
    oflline : Optional[bool]
        If True, offline reading is used with a local instance of DRUM, if
        availble. Default: False
    service_endpoint : Optional[str]
        Selects the TRIPS/DRUM web service endpoint to use. Is a choice between
        "drum" (default) and "drum-dev", a nightly build.

    Returns
    -------
    tp : TripsProcessor
        A TripsProcessor containing the extracted INDRA Statements
        in tp.statements.
    """
    if not offline:
        html = trips_client.send_query(text, service_endpoint)
        xml = trips_client.get_xml(html)
    else:
        if offline_reading:
            try:
                dr = DrumReader()
                if dr is None:
                    raise Exception('DrumReader could not be instantiated.')
            except BaseException as e:
                logger.error(e)
                logger.error('Make sure drum/bin/trips-drum is running in'
                              ' a separate process')
                return None
            try:
                dr.read_text(text)
                dr.start()
            except SystemExit:
                pass
            xml = dr.extractions[0]
        else:
            logger.error('Offline reading with TRIPS/DRUM not available.')
            logger.error('Error message was: %s' % offline_err)
            msg = """
                To install DRUM locally, follow instructions at
                https://github.com/wdebeaum/drum.
                Next, install the pykqml package either from pip or from
                https://github.com/bgyori/pykqml.
                Once installed, run drum/bin/trips-drum in a separate process.
                """
            logger.error(msg)
            return None
    if save_xml_name:
        trips_client.save_xml(xml, save_xml_name, save_xml_pretty)
    return process_xml(xml)


def process_xml(xml_string):
    """Return a TripsProcessor by processing a TRIPS EKB XML string.

    Parameters
    ----------
    xml_string : str
        A TRIPS extraction knowledge base (EKB) string to be processed.
        http://trips.ihmc.us/parser/api.html

    Returns
    -------
    tp : TripsProcessor
        A TripsProcessor containing the extracted INDRA Statements
        in tp.statements.
    """
    tp = TripsProcessor(xml_string)
    if tp.tree is None:
        return None
    tp.get_modifications_indirect()
    tp.get_activations_causal()
    tp.get_activations_stimulate()
    tp.get_complexes()
    tp.get_modifications()
    tp.get_active_forms()
    tp.get_active_forms_state()
    tp.get_activations()
    tp.get_translocation()
    tp.get_regulate_amounts()
    tp.get_degradations()
    tp.get_syntheses()
    tp.get_conversions()
    return tp
