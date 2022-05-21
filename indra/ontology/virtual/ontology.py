import requests
from ..ontology_graph import IndraOntology


class VirtualOntology(IndraOntology):
    """A virtual ontology class which uses a remote REST service to perform
    all operations. It is particularly useful if the host machine has limited
    resources and keeping the ontology graph in memory is not desirable.

    Parameters
    ----------
    url : str
        The base URL of the ontology graph web service.
    ontology : Optional[str]
        The identifier of the ontology recognized by the web service.
        Default: bio
    """
    def __init__(self, url, ontology='bio'):
        super().__init__()
        self.url = url
        self.ontology = ontology

    def initialize(self):
        self._initialized = True

    def child_rel(self, ns, id, rel_types):
        res = _send_request(self.url, 'child_rel',
                            ns=ns, id=id, rel_types=list(rel_types),
                            ontology=self.ontology)
        yield from (tuple(r) for r in res)

    def parent_rel(self, ns, id, rel_types):
        res = _send_request(self.url, 'parent_rel',
                            ns=ns, id=id, rel_types=list(rel_types),
                            ontology=self.ontology)
        yield from (tuple(r) for r in res)

    def get_node_property(self, ns, id, property):
        return _send_request(self.url, 'get_node_property',
                             ns=ns, id=id, property=property,
                             ontology=self.ontology)

    def get_id_from_name(self, ns, name):
        return _send_request(self.url, 'get_id_from_name',
                             ns=ns, name=name, ontology=self.ontology)


def _send_request(base_url, endpoint, **kwargs):
    url = '%s/%s' % (base_url, endpoint)
    res = requests.get(url, json=kwargs)
    return res.json()
