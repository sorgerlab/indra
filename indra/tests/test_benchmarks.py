from os.path import dirname, abspath, join
from indra.benchmarks import bioprocesses as bp
from indra.benchmarks import complexes as cp
from indra.benchmarks import phosphorylations as phos

eval_file = join(dirname(abspath(__file__)),
                 '../../models/assembly_eval/batch4/reach/' +
                 'reach_stmts_batch_4_eval.pkl')

def test_bioprocesses():
    """Smoke test to see if bioprocesses analysis works."""
    bp.analyze(eval_file)
