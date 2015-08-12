from __future__ import print_function, division

from sympy.logic.utilities import load_file
from sympy.logic import satisfiable
import time
import os
import sys

input_path = os.getcwd() + '/' + '/'.join(sys.argv[0].split('/')[:-1])

INPUT = [5 * i for i in range(2, 30)]
ALGORITHMS = ['dpll2', 'pycosat']
results = {}

for test in INPUT:
    results[test] = {}

for test in INPUT:
    for alg in ALGORITHMS:
        file_name = "%s/input/%d.cnf" % (input_path, test)
        theory = load_file(file_name)
        start = time.time()
        assert satisfiable(theory, algorithm=alg)
        end = time.time()
        results[test][alg] = end - start
        print("Test %d in time %.2f seconds for algorithm %s." %
            (test, end - start, alg))

print("problem," + ','.join(ALGORITHMS))

for test in INPUT:
    line = "%d" % test
    for alg in ALGORITHMS:
        line += ",%f" % results[test][alg]
    print(line)
