from __future__ import absolute_import, print_function, unicode_literals
import sys
from indra.benchmarks.bioprocesses import analyze

if __name__ == '__main__':
    # Get the command line args
    if len(sys.argv) < 2:
        print("Usage: %s stmts_file" % sys.argv[0])
        sys.exit()
    analyze(sys.argv[1])

