from indra.java_vm import autoclass, JavaException

class ReachReader(object):
    def __init__(self):
        self.api_ruler = None

    def get_api_ruler(self):
        if self.api_ruler is None:
            try:
                self.api_ruler =\
                    autoclass('edu.arizona.sista.reach.apis.ApiRuler')
            except JavaException:
                try:
                    autoclass('java.lang.String')
                except JavaException:
                    pass
                return None
        return self.api_ruler
