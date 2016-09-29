from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.tools.incremental_model import IncrementalModel
from indra.mechlinker import MechLinker
from indra.assemblers import EnglishAssembler


def print_linked_stmt(stmt):
    source_txts = []
    for source_stmt in stmt.source_stmts:
        source_txt = EnglishAssembler([source_stmt]).make_model()
        source_txts.append(source_txt)
    query_txt = EnglishAssembler([stmt.inferred_stmt]).make_model()
    final_txt =  'I know that '
    for i, t in enumerate(source_txts):
        final_txt += '(%d) %s ' % (i+1, t)
        if i < len(source_txts) -1:
            final_txt = final_txt[:-2] + ', and '
    final_txt += 'Is it therefore true that ' + query_txt[:-1] + '?'
    print(final_txt)
    return final_txt

if __name__ == '__main__':
    fname = 'models/rasmachine/rem/model.pkl'
    model = IncrementalModel(fname)
    model.preassemble()
    stmts = model.toplevel_stmts
    ml = MechLinker(stmts)
    linked_stmts = ml.link_statements()
    for stmt in linked_stmts:
        print_linked_stmt(stmt)
