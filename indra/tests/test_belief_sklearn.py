import random
from collections import defaultdict
from os.path import join, abspath, dirname
from indra.sources import signor
#from indra.belief.sklearn import SklearnBase
from indra.tools import assemble_corpus as ac

test_stmt_path = join(dirname(abspath(__file__)),
                      'belief_sklearn_test_stmts.pkl')
test_stmts = ac.load_statements(test_stmt_path)

# A set of statements derived from Signor used for testing purposes.
def _dump_test_data(filename, num_per_type=10):
    """Get corpus of statements for testing that has a range of stmt types."""
    sp = signor.process_from_web()
    # Group statements by type
    stmts_by_type = defaultdict(list)
    for stmt in sp.statements:
        stmts_by_type[stmt.__class__].append(stmt)
    # Sample statements of each type (without replacement)
    stmt_sample = []
    for stmt_type, stmt_list in stmts_by_type.items():
        if len(stmt_list) <= num_per_type:
            stmt_sample.extend(stmt_list)
        else:
            stmt_sample.extend(random.sample(stmt_list, num_per_type))
    ac.dump_statements(stmt_sample, filename)
    return stmt_sample



