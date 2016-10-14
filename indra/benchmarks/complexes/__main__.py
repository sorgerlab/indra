from __future__ import absolute_import, print_function, unicode_literals
import sys
from indra.benchmarks.complexes import analyze

if __name__ == '__main__':
    # Load the statements
    if len(sys.argv) < 2:
        print("Usage: %s reach_stmts_file" % sys.argv[0])
        sys.exit()
    results = analyze(sys.argv[1])
