from __future__ import print_function, division

from sympy.logic.utilities import load_file
from sympy.logic import satisfiable
import time
import os
import sys

input_path = os.path.dirname(__file__)

INPUT = [5 * i for i in range(2, 16)]
ALGORITHMS = ['dpll', 'dpll2']
results = {}

if __name__ == '__main__':
    for test in INPUT:
        results[test] = {}

    for test in INPUT:
        for alg in ALGORITHMS:
            file_name = os.path.join(input_path, 'input', '%s.cnf' % test)
            theory = load_file(file_name)
            start = time.time()
            if not satisfiable(theory, algorithm=alg):
                raise ValueError("Function returned false")
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
