from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import networkx
from indra.statements import *

class FigaroAssembler(object):
    def __init__(self, stmts, strengths=None, readout=None):
        self.statements = stmts
        self.strengths = strengths if strengths else {}
        self.BN = None
        self.magnitudes = {'small': 0.2, 'moderate': 0.5, 'huge': 2}
        self.readout = readout

    def make_model(self):
        self.BN = networkx.DiGraph()
        for stmt in self.statements:
            if isinstance(stmt, RegulateAmount):
                strength = self.strengths.get((stmt.subj.name, stmt.obj.name))
                if not strength:
                    coeff = 1
                else:
                    coeff = self.magnitudes[strength[1]] / \
                        self.magnitudes[strength[0]]
                if isinstance(stmt, DecreaseAmount):
                    pol = '-'
                elif isinstance(stmt, IncreaseAmount):
                    pol = '+'
                elif isinstance(stmt, Influence):
                    pol_stmt = stmt.overall_polarity()
                    pol = '+' if (pol_stmt is None or pol_stmt == '-1') \
                        else '-'
                else:
                    pol = '+'
                self.BN.add_edge(_n(stmt.subj.name), _n(stmt.obj.name),
                                 pol=pol, coeff=coeff)

    def print_model(self, fname=None):
        imports = ['language._', 'library.atomic._',
                   'library.compound._',
                   'library.atomic.continuous.Normal',
                   'algorithm.sampling.Importance']
        header = '\n'.join(['import com.cra.figaro.' + x for x in imports])
        wrapper = 'object IndraModel {\n' + \
            '%s\ndef main(args:Array[String]) = {%s}\n}'
        node_defs = []
        for node in self.BN.nodes():
            npar = self.BN.in_degree(node)
            if npar == 0:
                node_def = 'val %s = Normal(1,0.2)' % node
            elif npar == 1:
                node_def = 'val %s = Chain(%s, (v:Double) => Normal(v, 0.2))' % \
                    (node, list(self.BN.predecessors(node))[0])
            else:
                node_def_p = 'val %s = Chain(^^(%s), (v:(%s)) => Normal(%s, 0.2))'
                parents = ','.join(self.BN.predecessors(node))
                parent_types = ','.join(['Double'] * npar)
                terms = []
                for i, (p, _, meta) in enumerate(self.BN.in_edges(node, data=True)):
                    terms.append('%s %.2f * v._%d' % (meta['pol'], meta['coeff'], i+1))
                parent_fun = ' '.join(terms)
                node_def = node_def_p % (node, parents, parent_types, parent_fun)
            node_defs.append(node_def)
        if self.readout:
            main = 'var importance = Importance(100000, %s)\n' % self.readout
            main += 'importance.start()\n'
            main += 'val expv = importance.computeExpectation(%s, (v: Double) => v)\n' % \
                self.readout
            main += 'println(expv)'
        else:
            main = ''
        txt = header + '\n\n' + wrapper % ('\n'.join(node_defs), main)
        if fname:
            with open(fname, 'w') as fh:
                fh.write(txt)
        return txt

def _n(name):
    return name.replace(' ', '_')
