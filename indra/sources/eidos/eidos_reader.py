import os
import json

# Before the import, we have to deal with the CLASSPATH to avoid clashes
# with REACH.
def _set_classpath():
    clp = os.environ.get('CLASSPATH')
    eip = os.environ.get('EIDOSPATH')
    rep = os.environ.get('REACHPATH')
    clp_parts = clp.split(':') if clp else []
    new_clp_parts = []
    has_eidos = False
    # Look at all the parts of the CLASSPATH
    for part in clp_parts:
        # If REACH is on the CLASSPATH, remove it
        if not rep or os.path.abspath(part) != rep:
            new_clp_parts.append(part)
        # If Eidos is not on the CLASSPATH, add it
        if eip and os.path.abspath(part) == eip:
            has_eidos = True
    if eip and not has_eidos:
        new_clp_parts.append(eip)
    # Set the new CLASSPATH
    new_clp = ':'.join(new_clp_parts)
    os.environ['CLASSPATH'] = new_clp
_set_classpath()

from indra.java_vm import autoclass, JavaException


eidos_package = 'org.clulab.wm.eidos'


class EidosReader(object):
    """Reader object keeping an instance of the Eidos reader as a singleton.

    This allows the Eidos reader to need initialization when the first piece of
    text is read, the subsequent readings are done with the same
    instance of the reader and are therefore faster.

    Attributes
    ----------
    eidos_reader : org.clulab.wm.eidos.EidosSystem
        A Scala object, an instance of the Eidos reading system. It is
        instantiated only when first processing text.
    """

    def __init__(self):
        self.eidos_reader = None

    def process_text(self, text, format='json'):
        """Return a mentions JSON object given text.

        Parameters
        ----------
        text : str
            Text to be processed.
        format : str
            The format of the output to produce, one of "json" or "json_ld".
            Default: "json"

        Returns
        -------
        json_dict : dict
            A JSON object of mentions extracted from text.
        """
        if self.eidos_reader is None:
            eidos = autoclass(eidos_package + '.EidosSystem')
            self.eidos_reader = eidos(autoclass('java.lang.Object')())

        annot_doc = self.eidos_reader.extractFromText(text, False)
        if format == 'json':
            mentions = annot_doc.odinMentions()
            ser = autoclass(eidos_package + '.serialization.json.WMJSONSerializer')
            mentions_json = ser.toJsonStr(mentions)
        elif format == 'json_ld':
            # We need to get a Scala Seq of annot docs here
            ml = autoclass('scala.collection.mutable.MutableList')()
            ml.appendElem(annot_doc)
            jc = autoclass(eidos_package + '.serialization.json.JLDCorpus')
            corpus = jc(ml, None)
            mentions_json = corpus.toJsonStr()
        json_dict = json.loads(mentions_json)
        return json_dict

