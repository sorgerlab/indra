from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

if __name__ == '__main__':
    import pickle
    import sys

    with open(sys.argv[1]) as f:
        paper_stmts = pickle.load(f)

    for stmts in paper_stmts.values():
        for stmt in stmts:
            for ag in stmt.agent_list():
                if ag is None:
                    continue
                db_refs = ag.db_refs
                if ag.db_refs.get('HMDB'):
                    ag.db_refs['HMDB'] = 'HMDB%s' % ag.db_refs['HMDB']
                if ag.db_refs.get('CHEBI'):
                    ag.db_refs['CHEBI'] = 'CHEBI:%s' % ag.db_refs['CHEBI']
                if ag.db_refs.get('CHEMBL'):
                    ag.db_refs['CHEMBL'] = 'CHEMBL%s' % ag.db_refs['CHEMBL']
                if ag.db_refs.get('GO'):
                    ag.db_refs['GO'] = 'GO:%s' % ag.db_refs['GO']

    with open('fixed.pkl', 'w') as f:
        pickle.dump(paper_stmts, f)
