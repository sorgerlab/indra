from indra.statements import *


def process_pybel_graph(graph):
    proc = PybelProcessor(graph)
    proc.get_statements()
    return proc


class PybelProcessor(object):
    """Extract INDRA Statements from a PyBEL Graph.

    Parameters
    ----------
    graph : pybel.BELGraph
        PyBEL graph containing the BEL content.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of extracted INDRA Statements representing BEL Statements.
    """
    def __init__(self, graph):
        self.graph = graph


    def get_statements(self):
        for u, v, d in self.graph.edges_iter(data=True):
            pass
        st = Phosphorylation(Agent('A'), Agent('B'))
        self.statements = [st]
