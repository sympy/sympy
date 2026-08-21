"""
Execute like this:

$ python bin/generate_test_list.py
tests = [
    'sympy.concrete.tests',
    'sympy.core.tests',
    'sympy.functions.combinatorial.tests',
    'sympy.functions.elementary.tests',
    'sympy.functions.special.tests',
    'sympy.geometry.tests',
    'sympy.integrals.tests',
    'sympy.matrices.tests',
    'sympy.ntheory.tests',
    'sympy.numerics.tests',
    'sympy.parsing.tests',
    'sympy.physics.tests',
    'sympy.plotting.tests',
    'sympy.polynomials.tests',
    'sympy.printing.tests',
    'sympy.series.tests',
    'sympy.simplify.tests',
    'sympy.solvers.tests',
    'sympy.specfun.tests',
    'sympy.test_external',
    'sympy.utilities.tests',
    ]

"""

from __future__ import print_function

import os


def generate_test_list():
    """Generate the list of test package paths under the `sympy` tree.

    This uses os.walk so it works on Windows and POSIX. Any directory
    containing files named `test_*.py` is included as a dotted package
    (e.g. 'sympy.core.tests' or 'sympy').
    """
    g = []
    for root, dirs, files in os.walk('sympy'):
        for f in files:
            if f.startswith('test_') and f.endswith('.py'):
                pkg = root.replace(os.path.sep, '.')
                g.append(pkg)
                break

    g = list(set(g))
    g.sort()
    return g


if __name__ == '__main__':
    g = generate_test_list()
    print('tests = [')
    for x in g:
        print("    '%s'," % x)
    print(']')
