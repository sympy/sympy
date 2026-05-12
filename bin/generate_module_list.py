"""
Execute like this:

$ python bin/generate_module_list.py
modules = [
    'sympy.assumptions',
    'sympy.assumptions.handlers',
    'sympy.benchmarks',
    'sympy.calculus',
    'sympy.categories',
    'sympy.codegen',
    'sympy.combinatorics',
    'sympy.concrete',
    'sympy.core',
    'sympy.core.benchmarks',
    'sympy.crypto',
    'sympy.deprecated',
    'sympy.diffgeom',
    'sympy.external',
    'sympy.functions',
    'sympy.functions.combinatorial',
    'sympy.functions.elementary',
    'sympy.functions.elementary.benchmarks',
...
]

"""

from __future__ import print_function

import os


def generate_module_list():
    """Generate the list of package modules under the `sympy` tree.

    This implementation uses os.walk so it works consistently on Windows
    and POSIX systems. It returns a sorted list of dotted package paths
    (e.g. 'sympy.core', 'sympy.geometry'). The top-level 'sympy' package
    itself is excluded.
    """
    g = []
    for root, dirs, files in os.walk('sympy'):
        if '__init__.py' in files:
            # convert filesystem path to dotted module path
            pkg = root.replace(os.path.sep, '.')
            g.append(pkg)

    # remove the top-level package entry if present
    g = [i for i in g if not i.endswith('.tests')]
    if 'sympy' in g:
        g.remove('sympy')
    g = list(set(g))
    g.sort()
    return g


if __name__ == '__main__':
    g = generate_module_list()
    print('modules = [')
    for x in g:
        print("    '%s'," % x)
    print(']')
