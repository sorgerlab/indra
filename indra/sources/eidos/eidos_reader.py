import json
from indra.java_vm import autoclass, JavaException

class EidosReader(object):
    """Reader object keeping an instance of the Eidos reader as a singleton.

    This allows the Eidos reader to need initialization when the first piece of
    text is read, the subsequent readings are done with the same
    instance of the reader and are therefore faster.

    Attributes
    ----------
    eidos_reader : org.clulab.wm.AgroSystem
        A Scala object, an instance of the Eidos reading system. It is
        instantiated only when first processing text.
    """
    def __init__(self):
        self.eidos_reader = None

    def process_text(self, text):
        """Return a mentions JSON object given text.

        Parameters
        ----------
        text : str
            Text to be processed.

        Returns
        -------
        json_dict : dict
            A JSON object of mentions extracted from text.
        """
        if self.eidos_reader is None:
            eidos = autoclass('org.clulab.wm.AgroSystem')
            self.eidos_reader = eidos(autoclass('java.lang.Object')())

        mentions = self.eidos_reader.extractFrom(text)
        ser = autoclass('org.clulab.wm.serialization.json.WMJSONSerializer')
        mentions_json = ser.toJsonStr(mentions)
        json_dict = json.loads(mentions_json)
        return json_dict

