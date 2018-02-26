from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import logging

# Before the import, we have to deal with the CLASSPATH to avoid clashes
# with Eidos.
def _set_classpath():
    clp = os.environ.get('CLASSPATH')
    eip = os.environ.get('EIDOSPATH')
    rep = os.environ.get('REACHPATH')
    clp_parts = clp.split(':') if clp else []
    new_clp_parts = []
    has_reach = False
    # Look at all the parts of the CLASSPATH
    for part in clp_parts:
        # If Eidos is on the CLASSPATH, remove it
        if not eip or os.path.abspath(part) != eip:
            new_clp_parts.append(part)
        # If REACH is not on the CLASSPATH, add it
        if rep and os.path.abspath(part) == rep:
            has_reach = True
    if rep and not has_reach:
        new_clp_parts.append(rep)
    # Set the new CLASSPATH
    new_clp = ':'.join(new_clp_parts)
    os.environ['CLASSPATH'] = new_clp
_set_classpath()


from indra.java_vm import autoclass, JavaException

logger = logging.getLogger('reach_reader')

class ReachReader(object):
    """The ReachReader wraps a singleton instance of the REACH reader.

    This allows calling the reader many times without having to wait for it to
    start up each time.

    Attributes
    ----------
    api_ruler : org.clulab.reach.apis.ApiRuler
        An instance of the REACH ApiRuler class (java object).
    """
    def __init__(self):
        self.api_ruler = None

    def get_api_ruler(self):
        """Return the existing reader if it exists or launch a new one.

        Returns
        -------
        api_ruler : org.clulab.reach.apis.ApiRuler
            An instance of the REACH ApiRuler class (java object).
        """
        if self.api_ruler is None:
            try:
                self.api_ruler = \
                    autoclass('org.clulab.reach.export.apis.ApiRuler')
            except JavaException:
                # This second autoclass is needed because of a jnius
                # issue in which the first JavaException is not raised.
                try:
                    autoclass('java.lang.String')
                except JavaException as e:
                    logger.error(e)
                    pass
                return None
        return self.api_ruler
