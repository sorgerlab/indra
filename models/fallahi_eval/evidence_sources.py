from util import pklload
from collections import defaultdict
import indra.tools.assemble_corpus as ac


if __name__ == '__main__':
    # Load cached Statements just before going into the model
    stmts = pklload('pysb_stmts')

    # Start a dictionary for source counts
    sources_count = defaultdict(int)
    # Count statements according to sources of evidence
    for stmt in stmts:
        sources = tuple(sorted(list(set([ev.source_api for ev in stmt.evidence]))))
        sources_count[sources] += 1

    # Statements from databases only
    db_only = 0
    # Statements from reading only
    reading_only = 0
    # Statements from databases and reading
    mixture = 0
    # Database sources
    dbs = set(['bel', 'biopax', 'phosphosite', 'signor'])
    # Reader sources
    readers = set(['reach', 'trips', 'sparser', 'r3'])
    for k, v in sources_count.items():
        d = set(k).intersection(dbs)
        r = set(k).intersection(readers)
        if d and r:
            mixture += v
        if d and not r:
            db_only += v
        if r and not d:
            reading_only += v

    for k, v in sorted(sources_count.items(), key=lambda x: (len(x[0]), ','.join(sorted(x[0])))):
        sources_str = ','.join(k)
        line_str = sources_str + '\t' + str(v)
        print(line_str)
