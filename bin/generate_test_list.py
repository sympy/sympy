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

from glob import glob


def get_paths(level=15):
    """
    Generates a set of paths for testfiles searching.

    Example:
    >>> get_paths(2)
    ['sympy/test_*.py', 'sympy/*/test_*.py', 'sympy/*/*/test_*.py']
    >>> get_paths(6)
    ['sympy/test_*.py', 'sympy/*/test_*.py', 'sympy/*/*/test_*.py',
    'sympy/*/*/*/test_*.py', 'sympy/*/*/*/*/test_*.py',
    'sympy/*/*/*/*/*/test_*.py', 'sympy/*/*/*/*/*/*/test_*.py']

    """
    wildcards = ["/"]
    for i in range(level):
        wildcards.append(wildcards[-1] + "*/")
    p = ["sympy" + x + "test_*.py" for x in wildcards]
    return p

g = []
for x in get_paths():
    g.extend(glob(x))
g = [".".join(x.split("/")[:-1]) for x in g]
g = list(set(g))
g.sort()
print("tests = [")
for x in g:
    print("    '%s'," % x)
print("    ]")
